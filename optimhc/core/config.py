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
    "numProcess": 4,
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

    return config


class Config:
    """
    Configuration manager for optiMHC pipeline.
    Loads, validates, and provides access to configuration parameters from YAML or dict.
    """
    def __init__(self, config_source=None):
        """
        Initialize Config from a YAML file path or a dictionary.
        Args:
            config_source (str or dict, optional): Path to YAML file or dict with config.
        """
        if config_source is None:
            self._config = deepcopy(DEFAULT_CONFIG)
        elif isinstance(config_source, str):
            with open(config_source, "r") as f:
                user_config = yaml.safe_load(f)
            self._config = _deep_merge(DEFAULT_CONFIG, user_config)
        elif isinstance(config_source, dict):
            self._config = _deep_merge(DEFAULT_CONFIG, config_source)
        else:
            raise ValueError("Config source must be a file path, dict, or None.")
        # TODO: Add validation here

    def to_dict(self):
        return deepcopy(self._config)

    def save(self, path):
        with open(path, "w") as f:
            yaml.safe_dump(self._config, f)

    def __getitem__(self, key):
        return self._config[key]

    def __setitem__(self, key, value):
        self._config[key] = value

    def get(self, key, default=None):
        return self._config.get(key, default)

    def __contains__(self, key):
        return key in self._config

    def __repr__(self):
        return f"Config({self._config})"
