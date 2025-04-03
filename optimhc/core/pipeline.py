import os
import logging
import pandas as pd
from mokapot.model import PercolatorModel

from optimhc.parser import read_pin, read_pepxml
from optimhc.rescore import mokapot
from optimhc.rescore.model import XGBoostPercolatorModel, RandomForestPercolatorModel
from optimhc.visualization import (
    plot_feature_importance,
    visualize_target_decoy_features,
    visualize_feature_correlation,
    plot_qvalues
)
from optimhc.core.config import load_config
from optimhc.core.logging_helper import setup_loggers
from optimhc.core.feature_generation import import_generators, generate_features

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
        self.experiment = self.config.get('experimentName', 'optimhc_experiment') 
        self.output_dir = os.path.join(self.config['output_dir'], self.experiment)
        os.makedirs(self.output_dir, exist_ok=True)
        setup_loggers(os.path.join(self.output_dir, 'log'))

        
    def read_input(self):
        """
        Read input PSMs based on configuration.
        
        Returns:
            PsmContainer: An object containing the loaded PSMs.
        """
        input_type = self.config['inputType']
        input_files = self.config['inputFile']
        if not isinstance(input_files, list):
            input_files = [input_files]
        
        try:
            if input_type == 'pepxml':
                psms = read_pepxml(input_files, decoy_prefix=self.config['decoyPrefix'])
            elif input_type == 'pin':
                psms = read_pin(input_files)
            else:
                raise ValueError(f"Unsupported input type: {input_type}")
            return psms
        except Exception as e:
            logger.error(f"Failed to read input files: {e}")
            raise

            
    def rescore(self, psms, model_type=None, n_jobs=None, test_fdr=None, rescoring_features=None):
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
        rescore_config = self.config['rescore']
        
        # Use provided parameters or fall back to config values
        test_fdr = test_fdr if test_fdr is not None else rescore_config['testFDR']
        model_type = model_type if model_type is not None else rescore_config['model']
        n_jobs = n_jobs if n_jobs is not None else rescore_config['numJobs']
        
        if model_type == 'XGBoost':
            model = XGBoostPercolatorModel(n_jobs=n_jobs)
        elif model_type == 'RandomForest':
            model = RandomForestPercolatorModel(n_jobs=n_jobs)
        elif model_type == 'Percolator':
            model = PercolatorModel(n_jobs=n_jobs)
        else:
            model = PercolatorModel(n_jobs=n_jobs)
            
        kwargs = {}
        if rescoring_features is not None:
            kwargs['rescoring_features'] = rescoring_features
            
        results, models = mokapot.rescore(psms, model=model, test_fdr=test_fdr, **kwargs)
        return results, models
    
        
    def save_results(self, psms, results):
        """
        Save rescoring results and PSM data.
        
        Parameters:
            psms (PsmContainer): The PSM container object.
            results (mokapot.Results): The rescoring results.
        """
        results.to_txt(dest_dir=self.output_dir, file_root="optimhc", decoys=True)
        # psms.psms.to_csv(os.path.join(self.output_dir, 'psms.csv'), index=False)
        psms.write_pin(os.path.join(self.output_dir, f'{self.experiment}.optimhc.pin'))

        # Merge and save full information
        # df_combined = pd.concat([
        #     results.confidence_estimates['psms'],
        #     results.decoy_confidence_estimates['psms']
        # ], axis=0)

        # mokapot_columns = ['mokapot score', 'mokapot q-value', 'mokapot PEP']

        # df_full_information = pd.merge(
        #     df_combined[mokapot_columns + psms.identifier_columns],
        #     psms.psms,
        #     on=psms.identifier_columns,
        #     how='inner'
        # )

        # columns_order = [col for col in df_full_information.columns if col not in mokapot_columns] + mokapot_columns
        # df_full_information = df_full_information[columns_order]
        # df_full_information.to_csv(os.path.join(self.output_dir, 'results_psms.csv'), index=False)

        # logger.debug(f'Combined shape: {df_combined.shape}')
        # logger.debug(f'Full information shape: {df_full_information.shape}')

        
    def visualize_results(self, psms, results, models):
        """
        Visualize the results of the analysis pipeline.
        
        Parameters:
            psms (PsmContainer): The PSM container object.
            results (mokapot.Results): The rescoring results.
            models (list): The trained models used for rescoring.
        """
        if not self.config['visualization']:
            logger.info("Visualization is disabled. Skipping...")
            return
        
        fig_dir = os.path.join(self.output_dir, 'figures')

        os.makedirs(fig_dir, exist_ok=True)

        plot_qvalues(
            results,
            save_path=os.path.join(fig_dir, 'qvalues.png'),
            threshold=0.05,
        )

        plot_feature_importance(
            models,
            psms.rescoring_features,
            save_path=os.path.join(fig_dir, 'feature_importance.png')
        )
        visualize_target_decoy_features(
            psms,
            num_cols=4,
            save_path=os.path.join(fig_dir, 'target_decoy_histogram.png'),
        )
        visualize_feature_correlation(
            psms,
            save_path=os.path.join(fig_dir, 'feature_correlation.png'),
        )

        
    def run(self):
        """Run the complete pipeline."""
        logger.info("Starting analysis pipeline")
        
        # TODO: Dynamicaly import relevant generators (placeholder)
        # import_generators(self.config)
        
        # Read input PSMs
        psms = self.read_input()
        
        # Generate features
        generate_features(psms, self.config)
        
        # Rescore
        results, models = self.rescore(psms)
        
        # Save results
        self.save_results(psms, results)
        
        # Visualization
        self.visualize_results(psms, results, models)
        
        logger.info(f"Analysis pipeline completed, results saved to {self.output_dir}")

        
    def run_experiments(self):
        """Run experiments with different feature combinations."""
        logger.info("Starting experiment mode with multiple feature combinations")
        
        # Import relevant generators (placeholder)
        import_generators(self.config)
        
        # Read input PSMs
        psms = self.read_input()
        
        # Generate features
        generate_features(psms, self.config)
        
        # Save the initial PSMs data
        pin_path = os.path.join(self.output_dir, f'optimhc.{self.experiment}.pin')
        psms.write_pin(pin_path)
        # psms.psms.to_csv(os.path.join(self.output_dir, 'psms.csv'), index=False)
        
        # Run experiments
        experiment_configs = self.config.get('experiments', [])
        for i, exp_config in enumerate(experiment_configs):
            try:
                logger.info(f"Running experiment {i + 1}...")

                exp_name = exp_config.get('name', f'Experiment_{i + 1}')
                exp_dir = os.path.join(self.output_dir, exp_name)
                os.makedirs(exp_dir, exist_ok=True)

                source = exp_config.get('source', None)
                logger.info(f"Running experiment with sources: {source}")
                features = [feature for s in source for feature in psms.rescoring_features[s]]
                logger.info(f"Running experiment with features: {features}")
                
                n_jobs = self.config['global_parameters']['n_processes']
                logger.debug(f"n_jobs: {n_jobs}")
                
                try:
                    results, models = self.rescore(
                        psms, 
                        model_type=exp_config.get('model', None),
                        n_jobs=n_jobs,
                        test_fdr=0.01,
                        rescoring_features=features
                    )
                except Exception as e:
                    logger.error(f"Experiment {exp_name} failed during rescore: {e}")
                    continue  

                results.to_txt(dest_dir=exp_dir, decoys=True)

                for i, model in enumerate(models):
                    model.save(os.path.join(exp_dir, f'model.{i}'))
                    
                pin_path = os.path.join(exp_dir, f'{exp_name}.pin')
                psms.write_pin(pin_path, source=source)

                plot_qvalues(
                    results,
                    output_dir=exp_dir,
                    threshold=0.05
                )

                rescoring_features = {k: v for k, v in psms.rescoring_features.items() if k in source}
                plot_feature_importance(
                    models,
                    rescoring_features,
                    save_path=os.path.join(exp_dir, 'feature_importance.png')
                )
            except Exception as e:
                logger.error(f"Experiment {exp_name if 'exp_name' in locals() else i+1} encountered an unexpected error: {e}")
                continue
        
        logger.info("Experiments completed")
