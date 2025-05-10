import os
import sys
import logging
import click
import yaml
import json
from optimhc.core import Pipeline
from optimhc.core.config import Config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[logging.StreamHandler()]
)

logger = logging.getLogger(__name__)

@click.group()
@click.option('--logLevel', default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR)')
def cli(loglevel):
    """
    optiMHC: High-performance rescoring pipeline for immunopeptidomics data.

    Parameters
    ----------
    loglevel : str
        Logging level (DEBUG, INFO, WARNING, ERROR).
    """
    logging.getLogger().setLevel(getattr(logging, loglevel.upper(), logging.INFO))

def parse_cli_config(**kwargs):
    # Remove None values and build a config dict
    return {k: v for k, v in kwargs.items() if v is not None and v != ()}

@cli.command()
@click.option('--config', '-c', type=click.Path(exists=True), help='Path to YAML configuration file.')
@click.option('--inputType', type=click.Choice(['pepxml', 'pin']), help='Type of input file.')
@click.option('--inputFile', multiple=True, help='Path(s) to the input PSM file(s).')
@click.option('--decoyPrefix', help='Prefix used to identify decoy sequences.')
@click.option('--outputDir', help='Base directory for output files, logs, and figures.')
@click.option('--visualization/--no-visualization', default=None, help='Enable or disable visualization plots.')
@click.option('--removePreNxtAA/--keepPreNxtAA', default=None, help='Remove pre/post neighboring amino acids.')
@click.option('--numProcesses', type=int, help='Number of parallel processes to use.')
@click.option('--showProgress/--no-showProgress', default=None, help='Show progress information during execution.')
@click.option('--modificationMap', type=str, help='JSON string mapping modification masses to UNIMOD IDs.')
@click.option('--allele', multiple=True, help='Alleles for which predictions will be computed.')
@click.option('--featureGenerator', multiple=True, help='Feature generator configuration in JSON format.')
@click.option('--testFDR', type=float, help='FDR threshold for rescoring.')
@click.option('--model', type=str, help='Model to use for rescoring.')
@click.option('--numJobs', type=int, help='Number of parallel jobs for rescoring.')
def pipeline(**kwargs):
    """
    Run the standard optiMHC pipeline.

    Parameters
    ----------
    **kwargs
        Command-line options for pipeline configuration.
    """
    config_path = kwargs.pop('config', None)
    cli_config = parse_cli_config(**kwargs)

    # Parse JSON strings
    if 'modificationMap' in cli_config:
        cli_config['modificationMap'] = json.loads(cli_config['modificationMap'])
    if 'featureGenerator' in cli_config:
        fg = cli_config['featureGenerator']
        try:
            # Parse each feature generator config as JSON
            cli_config['featureGenerator'] = [json.loads(f) for f in fg]
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON format in feature generator config: {e}")
            raise
    if 'allele' in cli_config:
        cli_config['allele'] = list(cli_config['allele'])
    if 'inputFile' in cli_config:
        cli_config['inputFile'] = list(cli_config['inputFile'])

    if config_path:
        config = Config(config_path)
        # Override config with CLI config
        config._config = {**config._config, **cli_config}
    else:
        config = Config(cli_config)

    pipeline = Pipeline(config)
    pipeline.run()

@cli.command()
@click.option('--config', '-c', type=click.Path(exists=True), help='Path to YAML configuration file.')
@click.option('--input-type', type=click.Choice(['pepxml', 'pin']), help='Type of input file.')
@click.option('--input-file', multiple=True, help='Path(s) to the input PSM file(s).')
@click.option('--decoy-prefix', help='Prefix used to identify decoy sequences.')
@click.option('--output-dir', help='Base directory for output files, logs, and figures.')
@click.option('--visualization/--no-visualization', default=None, help='Enable or disable visualization plots.')
@click.option('--remove-pre-nxt-aa/--keep-pre-nxt-aa', default=None, help='Remove pre/post neighboring amino acids.')
@click.option('--num-processes', type=int, help='Number of parallel processes to use.')
@click.option('--show-progress/--no-show-progress', default=None, help='Show progress information during execution.')
@click.option('--modification-map', type=str, help='YAML/JSON string mapping modification masses to UNIMOD IDs.')
@click.option('--allele', multiple=True, help='Alleles for which predictions will be computed.')
@click.option('--feature-generator', multiple=True, help='Feature generator names (YAML/JSON for params).')
@click.option('--test-fdr', type=float, help='FDR threshold for rescoring.')
@click.option('--model', type=str, help='Model to use for rescoring.')
@click.option('--num-jobs', type=int, help='Number of parallel jobs for rescoring.')
def experiment(**kwargs):
    """
    Run optiMHC in experiment mode with multiple feature/model combinations.

    Parameters
    ----------
    **kwargs
        Command-line options for experiment configuration.
    """
    config_path = kwargs.pop('config', None)
    cli_config = parse_cli_config(**kwargs)
    import yaml
    if cli_config.get('modification_map'):
        cli_config['modificationMap'] = yaml.safe_load(cli_config.pop('modification_map'))
    if cli_config.get('feature_generator'):
        fg = cli_config.pop('feature_generator')
        try:
            cli_config['featureGenerator'] = yaml.safe_load('\n'.join(fg))
        except Exception:
            cli_config['featureGenerator'] = list(fg)
    if cli_config.get('allele'):
        cli_config['allele'] = list(cli_config['allele'])
    if cli_config.get('input_file'):
        cli_config['inputFile'] = list(cli_config['input_file'])
    if cli_config.get('output_dir'):
        cli_config['outputDir'] = cli_config.pop('output_dir')
    if cli_config.get('input_type'):
        cli_config['inputType'] = cli_config.pop('input_type')
    if cli_config.get('decoy_prefix'):
        cli_config['decoyPrefix'] = cli_config.pop('decoy_prefix')
    if cli_config.get('num_processes'):
        cli_config['numProcesses'] = cli_config.pop('num_processes')
    if cli_config.get('remove_pre_nxt_aa') is not None:
        cli_config['removePreNxtAA'] = cli_config.pop('remove_pre_nxt_aa')
    if cli_config.get('show_progress') is not None:
        cli_config['showProgress'] = cli_config.pop('show_progress')
    if cli_config.get('test_fdr'):
        cli_config.setdefault('rescore', {})['testFDR'] = cli_config.pop('test_fdr')
    if cli_config.get('model'):
        cli_config.setdefault('rescore', {})['model'] = cli_config.pop('model')
    if cli_config.get('num_jobs'):
        cli_config.setdefault('rescore', {})['numJobs'] = cli_config.pop('num_jobs')

    if config_path:
        config = Config(config_path)
        config._config = {**config._config, **cli_config}
    else:
        config = Config(cli_config)

    pipeline = Pipeline(config)
    pipeline.run_experiments()

if __name__ == "__main__":
    cli()
