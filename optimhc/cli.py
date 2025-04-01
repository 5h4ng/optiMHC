import os
import sys
import logging
import argparse

logger = logging.getLogger(__name__)

from optimhc.core import Pipeline


def main():
    """
    Command-line entry point for the analysis pipeline.
    Supports two modes:
    1. pipeline - Runs the standard pipeline
    2. experiment - Runs multiple experiments with different feature combinations
    """
    parser = argparse.ArgumentParser(
        description='Run Prophet analysis with the specified YAML configuration.'
    )
    parser.add_argument('config', type=str, help='Path to the YAML configuration file.')
    parser.add_argument('--mode', type=str, choices=['pipeline', 'experiment'], 
                        default='pipeline', help='Mode to run: standard pipeline or experiment mode')
    args = parser.parse_args()

    logger.info(f"Running analysis with configuration file: {args.config}, mode: {args.mode}")
    
    # Initialize the pipeline with the configuration file
    pipeline = Pipeline(args.config)
    
    # Run the pipeline in the selected mode
    if args.mode == 'experiment':
        logger.info("Starting experiment mode with multiple feature combinations")
        pipeline.run_experiments()
    else:
        logger.info("Starting standard pipeline mode")
        pipeline.run()

if __name__ == '__main__':
    main()
