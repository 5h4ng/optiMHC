# TODO: Validate the config file

import yaml
import logging
import os
from copy import deepcopy

logger = logging.getLogger(__name__)

# Default configuration with all possible parameters
DEFAULT_CONFIG = {
    'output_dir': './results',
    'log_file': None,  # Will be set to os.path.join(output_dir, 'log') if None
    'input_type': 'pepxml',
    'input_files': [],
    'decoy_prefix': 'DECOY_',
    'visualization': True,
    'score': None,
    'allele': [],
    'global_parameters': {
        'remove_modification': True,
        'remove_pre_nxt_aa': False,
        'n_processes': 1,
        'oxidation_tag': None,
        'show_progress': True
    },
    'rescore': {
        'test_fdr': 0.01,
        'model': 'percolator',
        'n_jobs': 1
    },
    'feature_generators': []
    # Additional feature generator default configs can be added as needed
}

def _deep_merge(default, override):
    """
    Deep merge two dictionaries. The override dictionary values take precedence.
    
    Parameters:
        default (dict): Default dictionary.
        override (dict): Dictionary with values to override defaults.
        
    Returns:
        dict: Merged dictionary.
    """
    result = deepcopy(default)
    
    if not isinstance(override, dict):
        return override
    
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
            
    return result

def load_config(config_path):
    """
    Load and parse a configuration file using YAML.
    Merges loaded config with default configuration.

    Parameters:
        config_path (str): Path to the YAML configuration file.

    Returns:
        dict: A dictionary containing all configurations.
    """
    logger.info(f"Loading configuration from {config_path}")
    with open(config_path, 'r') as f:
        user_config = yaml.safe_load(f)
    
    config = _deep_merge(DEFAULT_CONFIG, user_config)
    
    # Set log_file if not provided
    if config['log_file'] is None:
        config['log_file'] = os.path.join(config['output_dir'], 'log')
    
    # TODO: Add config validation here
    
    return config
