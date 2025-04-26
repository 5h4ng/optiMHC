# TODO: Validate the config file
# TODO: Create a config object to handle the config file

import yaml
import logging
import os
from copy import deepcopy

logger = logging.getLogger(__name__)

# Default configuration with all possible parameters
DEFAULT_CONFIG = {
    "outputDir": "./results",
    "inputType": "pepxml",
    "inputFile": [],
    "decoyPrefix": "DECOY_",
    "visualization": True,
    "saveModels": True,
    "allele": [],
    "numProcess": 32,
    "removePreNxtAA": False,
    "showProgress": True,
    "rescore": {"testFDR": 0.01, "model": "Percolator", "numJobs": 1},
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
    with open(config_path, "r") as f:
        user_config = yaml.safe_load(f)

    config = _deep_merge(DEFAULT_CONFIG, user_config)

    # TODO: Add config validation here

    return config
