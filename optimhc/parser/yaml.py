import logging
import yaml

logger = logging.getLogger(__name__)


def read_yaml(
    config_file: str
) -> dict:
    """
    Reads a YAML configuration file into a dictionary.

    Parameters:
        config_file (str): The file path to the YAML configuration file.
    
    Returns:
        dict: A dictionary containing the configuration data.
    """
    logger.info(f"Reading YAML configuration file: {config_file}")
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config