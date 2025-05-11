# psm_container.py

import logging
from typing import List, Optional, Union, Dict
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class PsmContainer:
    """
    A container for managing peptide-spectrum matches (PSMs) in immunopeptidomics rescoring pipelines.

    Parameters
    ----------
    psms : pd.DataFrame
        DataFrame containing the PSM data.
    label_column : str
        Column containing the label (True for target, False for decoy).
    scan_column : str
        Column containing the scan number.
    spectrum_column : str
        Column containing the spectrum identifier.
    ms_data_file_column : str
        Column containing the MS data file that the PSM originated from.
    peptide_column : str
        Column containing the peptide sequence.
    protein_column : str
        Column containing the protein accessions.
    rescoring_features : dict of str to list of str
        Dictionary of feature columns for rescoring.
    hit_rank_column : str, optional
        Column containing the hit rank.
    charge_column : str, optional
        Column containing the charge state.
    retention_time_column : str, optional
        Column containing the retention time.
    metadata_column : str, optional
        Column containing metadata.

    Attributes
    ----------
    psms : pd.DataFrame
        Copy of the DataFrame containing the PSM data.
    target_psms : pd.DataFrame
        DataFrame containing only target PSMs (label = True).
    decoy_psms : pd.DataFrame
        DataFrame containing only decoy PSMs (label = False).
    peptides : list of str
        List containing all peptides from the PSM data.
    columns : list of str
        List of column names in the PSM DataFrame.
    rescoring_features : dict of str to list of str
        Dictionary of rescoring feature columns in the PSM DataFrame.
    """

    def __init__(
        self,
        psms: pd.DataFrame,
        label_column: str,
        scan_column: str,
        spectrum_column: str,
        ms_data_file_column: str,
        peptide_column: str,
        protein_column: str,
        rescoring_features: Dict[str, List[str]],
        hit_rank_column: Optional[str] = None,
        charge_column: Optional[int] = None,
        retention_time_column: Optional[str] = None,
        metadata_column: Optional[str] = None,
    ):
        self._psms = psms.copy()
        self._psms.reset_index(drop=True, inplace=True)
        self.label_column = label_column
        self.scan_column = scan_column
        self.spectrum_column = spectrum_column
        self.ms_data_file_column = ms_data_file_column
        self.peptide_column = peptide_column
        self.protein_column = protein_column
        self.hit_rank_column = hit_rank_column
        self.retention_time_column = retention_time_column
        self.metadata_column = metadata_column
        self.rescoring_features = rescoring_features
        self.charge_column = charge_column
        # rescore result column
        self.rescore_result_column = None

        # check if the columns are in the dataframe
        def check_column(col):
            if col and col not in psms.columns:
                raise ValueError(f"Column '{col}' not found in PSM data.")

        check_column(label_column)
        check_column(scan_column)
        check_column(spectrum_column)
        check_column(ms_data_file_column)
        check_column(peptide_column)
        check_column(protein_column)
        check_column(hit_rank_column)
        check_column(retention_time_column)
        check_column(charge_column)

        if psms[label_column].nunique() == 1 and psms[label_column].iloc[0] == True:
            raise ValueError("All PSMs are labeled as target. No decoy PSMs found.")
        elif psms[label_column].nunique() == 1 and psms[label_column].iloc[0] == False:
            raise ValueError("All PSMs are labeled as decoy. No target PSMs found.")

        def check_metadata_column(col):
            # check the type is Dict[str, Dict[str, str]]
            if col and col not in psms.columns:
                raise ValueError(f"Column '{col}' not found in PSM data.")
            if not all(isinstance(x, dict) for x in self._psms[col]):
                raise ValueError(f"Column '{col}' must contain dictionaries.")

        if metadata_column:
            check_metadata_column(metadata_column)

        def check_rescoring_features(features: Dict[str, List[str]]):
            for key, cols in features.items():
                for col in cols:
                    if col not in psms.columns:
                        raise ValueError(
                            f"Column '{col}' not found in PSM data for feature '{key}'."
                        )

        check_rescoring_features(rescoring_features)
        
        # check if the number of decoy psms is not 0
        if len(self.decoy_psms) == 0:
            logger.error("No decoy PSMs found. Please check the decoy prefix.")
            raise ValueError("No decoy PSMs found.")

        logger.info("PsmContainer initialized with %d PSM entries.", len(self._psms))
        if self.ms_data_file_column:
            logger.info(
                "PSMs originated from %d MS data file(S).",
                len(self._psms[ms_data_file_column].unique()),
            )
        logger.info("target psms: %d", len(self.target_psms))
        logger.info("decoy psms: %d", len(self.decoy_psms))
        logger.info("unique peptides: %d", len(np.unique(self.peptides)))
        logger.info("rescoing features: %s", rescoring_features)
        

    @property
    def psms(self) -> pd.DataFrame:
        """
        Get a copy of the PSM DataFrame to prevent external modification.

        Returns
        -------
        pd.DataFrame
            A copy of the PSM DataFrame.
        """
        # TODO: return in a specific column order
        return self._psms.copy()

    def __len__(self) -> int:
        """
        Get the number of PSMs in the container.

        Returns
        -------
        int
            Number of PSMs.
        """
        return len(self._psms)

    @property
    def target_psms(self) -> pd.DataFrame:
        """
        Get a DataFrame containing only target PSMs.

        Returns
        -------
        pd.DataFrame
            DataFrame with only target PSMs (label = True).
        """
        return self._psms[self._psms[self.label_column] == True]

    @property
    def decoy_psms(self) -> pd.DataFrame:
        """
        Get a DataFrame containing only decoy PSMs.

        Returns
        -------
        pd.DataFrame
            DataFrame with only decoy PSMs (label = False).
        """
        return self._psms[self._psms[self.label_column] == False]

    @property
    def columns(self) -> List[str]:
        """
        Get the column names of the PSM DataFrame.

        Returns
        -------
        list of str
            List of column names.
        """
        return self._psms.columns

    @property
    def feature_columns(self) -> List[str]:
        """
        Get a list of all feature columns in the PSM DataFrame.

        Returns
        -------
        list of str
            List of feature column names.
        """
        return [col for cols in self.rescoring_features.values() for col in cols]

    @property
    def feature_sources(self) -> List[str]:
        """
        Get a list of all feature sources in the PSM DataFrame.

        Returns
        -------
        list of str
            List of feature source names.
        """
        return list(self.rescoring_features.keys())

    @property
    def peptides(self) -> List[str]:
        """
        Get the peptide sequences from the PSM data.

        Returns
        -------
        list of str
            List of peptide sequences.
        """
        return self._psms[self.peptide_column].tolist()

    @property
    def ms_data_files(self) -> List[str]:
        """
        Get the MS data files from the PSM data.

        Returns
        -------
        list of str
            List of MS data file names.
        """
        return self._psms[self.ms_data_file_column].tolist()

    @property
    def scan_ids(self) -> List[int]:
        """
        Get the scan numbers from the PSM data.

        Returns
        -------
        list of int
            List of scan numbers.
        """
        return self._psms[self.scan_column].tolist()

    @property
    def charges(self) -> List[int]:
        """
        Get the charge states from the PSM data.

        Returns
        -------
        list of int
            List of charge states.
        """
        return self._psms[self.charge_column].tolist()

    @property
    def metadata(self) -> pd.Series:
        """
        Get the metadata from the PSM data.

        Returns
        -------
        pd.Series
            Series containing metadata for each PSM.
        """
        return self._psms[self.metadata_column]

    @property
    def spectrum_ids(self) -> List[str]:
        """
        Get the spectrum identifiers from the PSM data.

        Returns
        -------
        list of str
            List of spectrum identifiers.
        """
        return self._psms[self.spectrum_column].tolist()

    @property
    def identifier_columns(self) -> List[str]:
        """
        Get the columns that uniquely identify each PSM.

        Returns
        -------
        list of str
            List of identifier column names.
        """
        return [
            self.scan_column,
            self.spectrum_column,
            self.peptide_column,
            self.protein_column,
        ]

    def copy(self) -> "PsmContainer":
        """
        Return a deep copy of the PsmContainer object.

        Returns
        -------
        PsmContainer
            A deep copy of the current PsmContainer.
        """
        import copy

        return copy.deepcopy(self)

    def __repr__(self) -> str:
        """
        Return a string representation of the PsmContainer.

        Returns
        -------
        str
            String summary of the PsmContainer.
        """
        return (
            f"PsmContainer with {len(self)} PSMs\n"
            f"\t - Target PSMs: {len(self.target_psms)}\n"
            f"\t - Decoy PSMs: {len(self.decoy_psms)}\n"
            f"\t - Unique Peptides: {len(np.unique(self.peptides))}\n"
            f"\t - Unique Spectra: {len(self._psms[self.spectrum_column].unique())}\n"
            f"\t - Rescoring Features: {self.rescoring_features}\n"
        )

    def drop_features(self, features: List[str]) -> None:
        """
        Drop specified features from the PSM DataFrame.

        Parameters
        ----------
        features : list of str
            List of feature column names to drop.

        Raises
        ------
        ValueError
            If any of the features do not exist in the DataFrame.
        """
        missing_features = [f for f in features if f not in self._psms.columns]
        if missing_features:
            raise ValueError(f"Features not found in PSM data: {missing_features}")

        self._psms.drop(columns=features, inplace=True)
        # Create a list of sources to update
        sources_to_update = []
        for source, cols in self.rescoring_features.items():
            self.rescoring_features[source] = [
                col for col in cols if col not in features
            ]
            if not self.rescoring_features[source]:
                sources_to_update.append(source)

        logger.info(
            f"Sources to be removed: {sources_to_update}. Because all the features are removed."
        )
        # Remove sources with no features left
        for source in sources_to_update:
            del self.rescoring_features[source]

    def drop_source(self, source: str) -> None:
        """
        Drop all features associated with a specific source from the PSM DataFrame.

        Parameters
        ----------
        source : str
            Name of the source to drop.

        Raises
        ------
        ValueError
            If the source does not exist in the rescoring features.
        """
        if source not in self.rescoring_features:
            raise ValueError(f"Source '{source}' not found in rescoring features.")
        self.drop_features(self.rescoring_features[source])

    def add_metadata(
        self,
        metadata_df: pd.DataFrame,
        psms_key: Union[str, List[str]],
        metadata_key: Union[str, List[str]],
        source,
    ) -> None:
        """
        Merge new metadata into the PSM DataFrame based on specified columns.
        Metadata from the specified source is stored as a nested dictionary inside the metadata column.

        Parameters
        ----------
        metadata_df : pd.DataFrame
            DataFrame containing new metadata to add.
        psms_key : str or list of str
            Column name(s) in the PSM data to merge on.
        metadata_key : str or list of str
            Column name(s) in the metadata data to merge on.
        source : str
            Name of the source of the new metadata.
        """
        if self.metadata_column is None:
            logger.info("No existing metadata column. Creating new metadata column.")
            self.metadata_column = "metadata"
            self._psms["metadata"] = [{} for _ in range(len(self._psms))]

        metadata_cols = [col for col in metadata_df.columns if col not in metadata_key]
        merged_df = self.psms.merge(
            metadata_df, left_on=psms_key, right_on=metadata_key, how="left"
        )
        if source in self._psms["metadata"]:
            logger.warning(f"{source} already exists in metadata. Overwriting.")
        for col in metadata_cols:
            merged_df["metadata"] = merged_df.apply(
                lambda row: {
                    **row["metadata"],
                    source: (
                        {col: row[col]}
                        if source not in row["metadata"]
                        else {**row["metadata"][source], col: row[col]}
                    ),
                },
                axis=1,
            )

        self._psms["metadata"] = merged_df["metadata"]

    def get_top_hits(self, n: int = 1):
        """
        Get the top n hits based on the hit rank column.
        If the hit rank column is not specified, returns the original PSMs.

        Parameters
        ----------
        n : int, optional
            The number of top hits to return. Default is 1.

        Returns
        -------
        PsmContainer
            A new PsmContainer object containing the top n hits.
        """
        if self.hit_rank_column is None:
            logger.warning("Rank column not specified. Return the original PSMs.")
            return self.copy()

        psms = self.copy()
        psms._psms = psms._psms[psms._psms[self.hit_rank_column] <= n]
        return psms

    def add_features(
        self,
        features_df: pd.DataFrame,
        psms_key: Union[str, List[str]],
        feature_key: Union[str, List[str]],
        source: str,
        suffix: Optional[str] = None,
    ) -> None:
        """
        Merge new features into the PSM DataFrame based on specified columns.

        Parameters
        ----------
        features_df : pd.DataFrame
            DataFrame containing new features to add.
        psms_key : str or list of str
            Column name(s) in the PSM data to merge on.
        feature_key : str or list of str
            Column name(s) in the features data to merge on.
        source : str
            Name of the source of the new features.
        suffix : str, optional
            Suffix to add to the new columns if there's a name conflict.
        """
        if isinstance(psms_key, str):
            psms_key = [psms_key]

        if isinstance(feature_key, str):
            feature_key = [feature_key]

        new_feature_cols = [
            col for col in features_df.columns if col not in feature_key
        ]

        for cols in new_feature_cols:
            if cols in self._psms.columns:
                logger.warning(f"Column '{cols}' already exists in PSM data.")
                if suffix is None:
                    logger.warning("No suffix provided. Using default suffix ")
                    raise ValueError("Duplicate columns exist. No suffix provided.")
                else:
                    logger.warning(
                        f"Suffix '{suffix}' provided. Using suffix '{suffix}'."
                    )
        logger.info(f"Adding {len(new_feature_cols)} new features from {source}.")

        if not new_feature_cols:
            logger.warning("No new features to add.")
            raise ValueError("No new features to add.")

        if source in self.rescoring_features:
            logger.warning(
                f"{source} already exists in rescoring features. Overwriting."
            )
            self.drop_source(source)

        # TODO: reluctant logic
        if suffix is None:
            suffixes = ("", "")
        else:
            suffixes = ("", suffix)

        self.rescoring_features[source] = [
            col + suffixes[1] for col in new_feature_cols
        ]
        features_df = features_df.rename(
            columns={col: col + suffixes[1] for col in new_feature_cols}
        )
        original_len = len(self._psms)
        # avoid merge the right key to the psms
        self._psms = self._psms.merge(
            features_df, left_on=psms_key, right_on=feature_key, how="left"
        )

        if feature_key != psms_key:
            cols_to_drop = [
                col
                for col in feature_key
                if col not in psms_key and col in self._psms.columns
            ]
            if cols_to_drop:
                logger.debug(
                    f"Dropping columns from feature_key not in psms_key: {cols_to_drop}"
                )
                self._psms.drop(columns=cols_to_drop, inplace=True)

        if len(self._psms) != original_len:
            raise ValueError(
                "Merging features resulted in a change in the number of PSMs. Check for duplicate keys."
            )

    def add_features_by_index(
        self, features_df: pd.DataFrame, source: str, suffix: Optional[str] = None
    ) -> None:
        """
        Merge new features into the PSM DataFrame based on the DataFrame index.

        Parameters
        ----------
        features_df : pd.DataFrame
            DataFrame containing new features to add.
        source : str
            Name of the source of the new features.
        suffix : str, optional
            Suffix to add to the new columns if there's a name conflict.
        """
        new_feature_cols = [col for col in features_df.columns]
        for col in new_feature_cols:
            if col in self._psms.columns:
                logger.warning(f"Column '{col}' already exists in PSM data.")
                if suffix is None:
                    logger.warning("No suffix provided. Using default suffix.")
                    raise ValueError("Duplicate columns exist. No suffix provided.")
                else:
                    logger.warning(
                        f"Suffix '{suffix}' provided. Using suffix '{suffix}'."
                    )

        logger.info(
            f"Adding {len(new_feature_cols)} new features from {source} by index."
        )

        if not new_feature_cols:
            logger.warning("No new features to add.")
            raise ValueError("No new features to add.")

        if source in self.rescoring_features:
            logger.warning(
                f"{source} already exists in rescoring features. Overwriting."
            )
            self.drop_source(source)

        if suffix is None:
            suffixes = ("", "")
        else:
            suffixes = ("", suffix)

        self.rescoring_features[source] = [
            col + suffixes[1] for col in new_feature_cols
        ]
        features_df.rename(
            columns={col: col + suffixes[1] for col in new_feature_cols}, inplace=True
        )
        original_len = len(self._psms)
        self._psms = self._psms.merge(
            features_df,
            left_index=True,
            right_index=True,
            how="left",  # Perform a left join to preserve all original PSM data
        )

        # Ensure that the merge did not change the number of rows in the PSM DataFrame
        if len(self._psms) != original_len:
            raise ValueError(
                "Merging features resulted in a change in the number of PSMs. Check for duplicate indices."
            )

    def add_results(
        self,
        results_df: pd.DataFrame,
        psms_key: Union[str, List[str]],
        result_key: Union[str, List[str]],
    ) -> None:
        """
        Add results of rescore engine to the PSM DataFrame based on specified columns.

        Parameters
        ----------
        results_df : pd.DataFrame
            DataFrame containing new results to add.
        psms_key : str or list of str
            Column name(s) in the PSM data to merge on.
        result_key : str or list of str
            Column name(s) in the results data to merge on.
        """
        if self.rescore_result_column is not None:
            logger.warning("Rescore result column already exists. Overwriting.")

        if set(self._psms.columns) & set(results_df.columns):
            raise ValueError(
                "Duplicate columns exist. Please rename the columns in the results data."
            )

        self.rescore_result_column = result_key
        self._psms = self._psms.merge(
            results_df,
            left_on=psms_key,
            right_on=result_key,
            how="left",
            validate="one_to_one",
        )
        self._psms.drop(columns=result_key, inplace=True)
        logger.info(f"Added rescore results to PSM data.")

    def write_pin(self, output_file: str, source: List[str] = None) -> None:
        """
        Write the PSM data to a Percolator input (PIN) file.

        Percolator accepts input in a simple tab-delimited format where each row contains features associated with a single PSM:
            id <tab> label <tab> scannr <tab> feature1 <tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM
        With header:
            SpecID <tab> Label <tab> ScanNr <tab> Feature1 <tab> ... <tab> FeatureN <tab> Peptide <tab> Proteins

        Parameters
        ----------
        output_file : str
            The path to the output PIN file.
        source : list of str, optional
            List of feature sources to include. If None, include all sources.

        Returns
        -------
        pd.DataFrame
            The DataFrame written to the PIN file.
        """
        df = self._psms.copy()
        # Check if the label column is str
        # Case1: label column is str
        if df[self.label_column].dtype == "str":
            df["PercolatorLabel"] = df[self.label_column].map({"True": 1, "False": -1})
        # Case2: label column is bool
        elif df[self.label_column].dtype == "bool":
            df["PercolatorLabel"] = df[self.label_column].map({True: 1, False: -1})
        else:
            # try to convert to bool
            logger.warning("Label column is not str or bool. Converting to bool.")
            df["PercolatorLabel"] = (
                df[self.label_column].astype(bool).map({True: 1, False: -1})
            )
        logger.info("Writing PIN file to %s", output_file)

        feature_cols = []
        if source is None:
            for _, cols in self.rescoring_features.items():
                feature_cols.extend(cols)
        else:
            for s in source:
                if s not in self.rescoring_features:
                    raise ValueError(f"Source '{s}' not found in rescoring features.")
                feature_cols.extend(self.rescoring_features[s])

        pin_df = pd.DataFrame()
        pin_df["SpecID"] = df[self.spectrum_column]
        pin_df["Label"] = df["PercolatorLabel"]
        pin_df["ScanNr"] = df[self.scan_column]
        for col in feature_cols:
            pin_df[col] = df[col]

        pin_df["Peptide"] = df[self.peptide_column]
        pin_df["Proteins"] = df[self.protein_column].apply(
            lambda x: ";".join(x) if isinstance(x, (list, tuple)) else x
        )
        pin_df.to_csv(output_file, sep="\t", index=False)
        logger.info("PIN file written to %s", output_file)

        return pin_df
