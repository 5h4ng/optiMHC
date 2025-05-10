# __init__.py for optimhc.core
# NOTE: To speed up CLI help and import, avoid importing heavy dependencies (deeplc, mhcflurry, etc.) at the top level.
# Import them only inside the functions/classes that require them.

from optimhc.core.pipeline import Pipeline

__all__ = ["Pipeline"]
