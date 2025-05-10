# rescore/mokapot.py

import logging
from typing import List, Dict
import mokapot
from mokapot import LinearPsmDataset
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


def rescore(
    psms: PsmContainer,
    model=None,
    rescoring_features: List[str] = None,
    test_fdr: float = 0.01,
    **kwargs,
):
    """
    Rescore PSMs using mokapot.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing PSM data.
    model : object, optional
        A trained model for rescoring PSMs.
    rescoring_features : List[str], optional
        A list of feature names to use for rescoring.
    test_fdr : float, optional
        The FDR threshold for testing the model. Default is 0.01.
    **kwargs : dict
        Additional keyword arguments for mokapot.brew.

    Returns
    -------
    tuple
        A tuple containing:
        - Confidence object or list of Confidence objects:
          An object or a list of objects containing the confidence estimates at various levels
          (i.e. PSMs, peptides) when assessed using the learned score. If a list, they will be
          in the same order as provided in the psms parameter.
        - list of Model objects:
          The learned Model objects, one for each fold.

    Notes
    -----
    This function:
    1. Converts the PsmContainer to a mokapot dataset
    2. Runs mokapot.brew with the specified parameters
    3. Returns the results and models
    """
    psms = convert_to_mokapot_dataset(psms, rescoring_features=rescoring_features)
    logger.info("Rescoring PSMs with mokapot.")
    results, models = mokapot.brew(psms, model=model, test_fdr=test_fdr, **kwargs)
    return results, models


def convert_to_mokapot_dataset(
    psms: PsmContainer, rescoring_features: List[str] = None
) -> LinearPsmDataset:
    """
    Convert a PsmContainer to a LinearPsmDataset for use with mokapot.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing PSM data.
    rescoring_features : List[str], optional
        A list of feature names to use for rescoring.
        If not provided, uses all features from the PsmContainer.

    Returns
    -------
    LinearPsmDataset
        A LinearPsmDataset object for use with mokapot.

    Raises
    ------
    ValueError
        If any of the specified rescoring features are not found in the PSM data.

    Notes
    -----
    This function:
    1. Extracts all features from the PsmContainer
    2. Validates the specified rescoring features
    3. Creates a LinearPsmDataset with the appropriate columns and data
    """

    # rescoring_features: Dict[str, List[str]] -> Tuple[str]
    feature_columns = [
        col for features in psms.rescoring_features.values() for col in features
    ]

    if rescoring_features is None:
        rescoring_features = feature_columns
    else:
        for feature in rescoring_features:
            if feature not in feature_columns:
                raise ValueError(f"Feature '{feature}' not found in the PSM data.")

    dataset = LinearPsmDataset(
        psms.psms,
        target_column=psms.label_column,
        spectrum_columns=psms.spectrum_column,
        peptide_column=psms.peptide_column,
        protein_column=psms.protein_column,
        feature_columns=rescoring_features,
        scan_column=psms.scan_column,
        rt_column=psms.retention_time_column,
    )

    return dataset
