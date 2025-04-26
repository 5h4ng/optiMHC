import os
import logging
import gc
import pandas as pd
from multiprocessing import Process
from mokapot.model import PercolatorModel

from optimhc.parser import read_pin, read_pepxml
from optimhc.rescore import mokapot
from optimhc.rescore.model import XGBoostPercolatorModel, RandomForestPercolatorModel
from optimhc.visualization import (
    plot_feature_importance,
    visualize_target_decoy_features,
    visualize_feature_correlation,
    plot_qvalues,
)
from optimhc.core.config import load_config
from optimhc.core.logging_helper import setup_loggers
from optimhc.core.feature_generation import generate_features

logger = logging.getLogger(__name__)


class Pipeline:
    """Pipeline class that encapsulates the entire data processing workflow."""

    def __init__(self, config_path):
        """
        Initialize the pipeline with a configuration file.

        Parameters:
            config_path (str): Path to the YAML configuration file.
        """
        self.config = load_config(config_path)
        self.experiment = self.config.get("experimentName", "optimhc_experiment")
        self.output_dir = os.path.join(self.config["outputDir"], self.experiment)
        os.makedirs(self.output_dir, exist_ok=True)
        setup_loggers(os.path.join(self.output_dir, "log"))

        self.visualization_enabled = self.config.get("visualization", True)
        self.save_models = self.config.get("saveModels", True)
        self.test_fdr = self.config.get("rescore", {}).get("testFDR", 0.01)
        self.model_type = self.config.get("rescore", {}).get("model", "Percolator")
        self.n_jobs = self.config.get("rescore", {}).get("numJobs", 1)

    def read_input(self):
        """
        Read input PSMs based on configuration.

        Returns:
            PsmContainer: An object containing the loaded PSMs.
        """
        input_type = self.config["inputType"]
        input_files = self.config["inputFile"]
        if not isinstance(input_files, list):
            input_files = [input_files]

        try:
            if input_type == "pepxml":
                psms = read_pepxml(input_files, decoy_prefix=self.config["decoyPrefix"])
            elif input_type == "pin":
                psms = read_pin(input_files)
            else:
                raise ValueError(f"Unsupported input type: {input_type}")
            return psms
        except Exception as e:
            logger.error(f"Failed to read input files: {e}")
            raise

    def _generate_features(self, psms):
        """
        Generate features for PSMs.

        Parameters:
            psms (PsmContainer): The PSM container object.

        Returns:
            PsmContainer: The PSM container with generated features.
        """
        generate_features(psms, self.config)
        return psms

    def rescore(
        self, psms, model_type=None, n_jobs=None, test_fdr=None, rescoring_features=None
    ):
        """
        Perform rescoring on the PSMs.

        Parameters:
            psms (PsmContainer): The PSM container object.
            model_type (str, optional): Type of model to use ('xgboost', 'random_forest', 'percolator').
            n_jobs (int, optional): Number of parallel jobs.
            test_fdr (float, optional): FDR threshold for testing.
            rescoring_features (list, optional): Features to use for rescoring.

        Returns:
            tuple: (results, models) - The rescoring results and trained models.
        """
        # Use provided parameters or fall back to config values
        test_fdr = test_fdr if test_fdr is not None else self.test_fdr
        model_type = model_type if model_type is not None else self.model_type
        n_jobs = n_jobs if n_jobs is not None else self.n_jobs

        if model_type == "XGBoost":
            model = XGBoostPercolatorModel(n_jobs=n_jobs)
        elif model_type == "RandomForest":
            model = RandomForestPercolatorModel(n_jobs=n_jobs)
        elif model_type == "Percolator":
            model = PercolatorModel(n_jobs=n_jobs)
        else:
            model = PercolatorModel(n_jobs=n_jobs)

        kwargs = {}
        if rescoring_features is not None:
            kwargs["rescoring_features"] = rescoring_features

        results, models = mokapot.rescore(
            psms, model=model, test_fdr=test_fdr, **kwargs
        )
        return results, models

    def save_results(self, psms, results, models, output_dir=None, file_root="optimhc"):
        """
        Save rescoring results and PSM data.

        Parameters:
            psms (PsmContainer): The PSM container object.
            results (mokapot.Results): The rescoring results.
            models (list): The trained models.
            output_dir (str, optional): Directory to save results to (defaults to self.output_dir).
            file_root (str, optional): Root name for output files.
        """
        output_dir = output_dir if output_dir is not None else self.output_dir

        results.to_txt(dest_dir=output_dir, file_root=file_root, decoys=True)
        psms.write_pin(os.path.join(output_dir, f"{file_root}.pin"))

        if self.save_models:
            model_dir = os.path.join(output_dir, "models")
            os.makedirs(model_dir, exist_ok=True)
            logger.info(f"Saving models to {model_dir}")
            for i, model in enumerate(models):
                model.save(os.path.join(model_dir, f"{file_root}.model{i}"))

    def visualize_results(self, psms, results, models, output_dir=None, sources=None):
        """
        Visualize the results of the analysis pipeline.

        Parameters:
            psms (PsmContainer): The PSM container object.
            results (mokapot.Results): The rescoring results.
            models (list): The trained models used for rescoring.
            output_dir (str, optional): Directory to save visualizations to.
            sources (list, optional): Sources of features to include in visualizations.
        """
        if not self.visualization_enabled:
            logger.info("Visualization is disabled. Skipping...")
            return

        output_dir = output_dir if output_dir is not None else self.output_dir
        fig_dir = os.path.join(output_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        plot_qvalues(
            results,
            save_path=os.path.join(fig_dir, "qvalues.png"),
            threshold=0.05,
        )

        if sources:
            rescoring_features = {
                k: v for k, v in psms.rescoring_features.items() if k in sources
            }
        else:
            rescoring_features = psms.rescoring_features

        plot_feature_importance(
            models,
            rescoring_features,
            save_path=os.path.join(fig_dir, "feature_importance.png"),
        )
        visualize_feature_correlation(
            psms,
            save_path=os.path.join(fig_dir, "feature_correlation.png"),
        )
        visualize_target_decoy_features(
            psms,
            num_cols=4,
            save_path=os.path.join(fig_dir, "target_decoy_histogram.png"),
        )

    def _run_single_experiment(self, psms, exp_config, exp_name, exp_dir):
        """
        Run a single experiment with the specified configuration.

        Parameters:
            psms (PsmContainer): The PSM container object.
            exp_config (dict): Configuration for this experiment.
            exp_name (str): Name of the experiment.
            exp_dir (str): Directory to save experiment results.

        Returns:
            bool: True if experiment succeeded, False otherwise.
        """
        try:
            os.makedirs(exp_dir, exist_ok=True)

            source = exp_config.get("source", None)
            model_type = exp_config.get("model", self.model_type)
            n_jobs = exp_config.get("numJobs", self.n_jobs)

            logger.info(f"Running experiment '{exp_name}' with sources: {source}")

            # Generate list of features based on the provided sources
            features = []
            if source:
                for s in source:
                    features.extend(psms.rescoring_features.get(s, []))
            logger.info(f"Features used in experiment '{exp_name}': {features}")

            results, models = self.rescore(
                psms,
                model_type=model_type,
                n_jobs=n_jobs,
                test_fdr=self.test_fdr,
                rescoring_features=features,
            )

            self.save_results(
                psms, results, models, output_dir=exp_dir, file_root=exp_name
            )

            fig_dir = os.path.join(exp_dir, "figures")

            plot_qvalues(
                results,
                save_path=os.path.join(fig_dir, "qvalues.png"),
                threshold=0.05,
            )

            plot_feature_importance(
                models,
                rescoring_features={
                    k: v for k, v in psms.rescoring_features.items() if k in source
                },
                save_path=os.path.join(fig_dir, "feature_importance.png"),
            )

            return True

        except Exception as e:
            logger.error(f"Experiment '{exp_name}' failed: {e}")
            return False

        finally:
            # Explicit resource release to free up memory after each experiment
            try:
                del results
                del models
            except Exception:
                pass
            gc.collect()

    def run(self):
        """Run the complete pipeline."""
        logger.info("Starting analysis pipeline")

        psms = self.read_input()
        psms = self._generate_features(psms)
        results, models = self.rescore(psms)
        self.save_results(psms, results, models)
        self.visualize_results(psms, results, models)

        logger.info(f"Analysis pipeline completed, results saved to {self.output_dir}")
        return psms, results, models

    def run_experiments(self):
        """
        Run experiments with different feature combinations using multiprocessing.

        Each experiment is executed in its own process for complete resource isolation.

        Attention: The config must
        """
        logger.info("Starting experiment mode with multiple feature combinations")

        psms = self.read_input()
        psms = self._generate_features(psms)

        # Save the generated pin file for reference
        pin_path = os.path.join(self.output_dir, f"optimhc.{self.experiment}.pin")
        psms.write_pin(pin_path)
        fig_summary_dir = os.path.join(self.output_dir, "figures")
        os.makedirs(fig_summary_dir, exist_ok=True)
        visualize_feature_correlation(
            psms,
            save_path=os.path.join(fig_summary_dir, "feature_correlation.png"),
        )
        # visualize_target_decoy_features(
        #     psms,
        #     num_cols=4,
        #     save_path=os.path.join(fig_summary_dir, 'target_decoy_histogram.png'),
        # )

        experiment_configs = self.config.get("experiments", [])
        processes = []
        for i, exp_config in enumerate(experiment_configs):
            exp_name = exp_config.get("name", f"Experiment_{i + 1}")
            exp_dir = os.path.join(self.output_dir, exp_name)

            logger.info(f"Starting experiment '{exp_name}' in a separate process")
            p = Process(
                target=self._run_single_experiment,
                args=(psms, exp_config, exp_name, exp_dir),
            )
            p.start()
            processes.append(p)

        # Wait for all experiment processes to finish
        for p in processes:
            p.join()

        logger.info("All experiments completed")
