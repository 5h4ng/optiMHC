import logging

def setup_loggers(log_file=None):
    """
    Create or update all loggers so that each logger has a StreamHandler and optionally a FileHandler.
    This ensures all log messages are displayed in the console and optionally saved to a file.

    Parameters:
        log_file (str): Path to the log file. If None, no file logging is set up.
    """
    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for lg in loggers:
        lg.disabled = False
        has_stream_handler = any(isinstance(handler, logging.StreamHandler) for handler in lg.handlers)
        if not has_stream_handler:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            lg.addHandler(console_handler)

        if log_file:
            has_file_handler = any(isinstance(handler, logging.FileHandler) for handler in lg.handlers)
            if not has_file_handler:
                file_handler = logging.FileHandler(log_file, mode='a')
                file_handler.setLevel(logging.INFO)
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                file_handler.setFormatter(formatter)
                lg.addHandler(file_handler)

        lg.propagate = False

        if lg.name.startswith('optimhc'):
            lg.disabled = False

    root_logger = logging.getLogger()
    root_logger.disabled = False

def debug_logging():
    """
    Print debugging information for all loggers that start with 'optimhc' and
    the root logger. This helps verify that logger configurations are set properly.
    """
    print("\n=== Debugging Loggers ===\n")
    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict.keys()]
    for lg in loggers:
        if lg.name.startswith("optimhc"):
            print(f"Logger Name: {lg.name}")
            print(f"  - Effective Level: {logging.getLevelName(lg.getEffectiveLevel())}")
            print(f"  - Explicit Level: {logging.getLevelName(lg.level)} (default: NOTSET)")
            print(f"  - Propagate: {lg.propagate}")
            print(f"  - Disabled: {lg.disabled}")

            if lg.handlers:
                for handler in lg.handlers:
                    print(f"    Handler: {type(handler).__name__}")
                    print(f"      - Level: {logging.getLevelName(handler.level)}")
                    print(f"      - Formatter: {handler.formatter}")
                    if isinstance(handler, logging.FileHandler):
                        print(f"      - Log File: {handler.baseFilename}")
                    print(f"      - Stream: {getattr(handler, 'stream', None)}")
            else:
                print(f"    No handlers attached to the logger.")
            print("")

    root_logger = logging.getLogger()
    print(f"Root Logger:")
    print(f"  - Level: {logging.getLevelName(root_logger.level)}")
    print(f"  - Handlers: {len(root_logger.handlers)}")
    for handler in root_logger.handlers:
        print(f"    Handler: {type(handler).__name__}")
        print(f"      - Level: {logging.getLevelName(handler.level)}")
        print(f"      - Formatter: {handler.formatter}")
        if isinstance(handler, logging.FileHandler):
            print(f"      - Log File: {handler.baseFilename}")
        print(f"      - Stream: {getattr(handler, 'stream', None)}")
    print("\n=== End of Logger Debugging ===\n")
