# feature_generator/PWM.py

# TODO: extract best allele for each peptide
import os
import pandas as pd
import numpy as np
from typing import Optional, Dict, Union, List, Tuple
from optimhc.feature_generator.base_feature_generator import BaseFeatureGenerator
import logging
from optimhc import utils

logger = logging.getLogger(__name__)

# n_flank_pwm and c_flank_pwm are used for MHC class II PWM calculation.
# The original PPM is recorded in data/PWMs/II/n_flank_ppm and data/PWMs/II/c_flank_ppm,
# which is converted to PWM by taking the background frequency as 0.05.

c_flank_pwm_data = {
    "A": [0.891194, 0.599125, 0.567353],
    "C": [-3.952777, -4.345220, -4.612584],
    "D": [0.291173, 0.550133, 0.325690],
    "E": [0.687212, 0.662834, 0.834717],
    "F": [-1.250652, -0.784627, -1.139232],
    "G": [0.509354, -0.368370, 0.919885],
    "H": [-0.808229, -0.836628, -0.591508],
    "I": [-0.046196, -0.034123, -1.564360],
    "K": [0.452471, 0.665617, 1.086677],
    "L": [0.522681, 0.382607, 0.208291],
    "M": [-2.131278, -2.144693, -2.008653],
    "N": [-0.208044, -0.699085, 0.022668],
    "P": [0.417673, 1.179052, -0.018850],
    "Q": [-0.033656, -0.051822, 0.274208],
    "R": [0.535173, 0.397917, 0.924583],
    "S": [0.542669, 0.255087, 0.360228],
    "T": [-0.359617, -0.322741, -0.028565],
    "V": [0.507861, 0.582748, -0.584754],
    "W": [-3.432776, -3.024952, -3.391619],
    "Y": [-1.144700, -0.725718, -1.221789],
    # Add a pseudo-entry for 'X' with 0 contribution
    "X": [0.0, 0.0, 0.0],
}

c_flank_pwm = pd.DataFrame.from_dict(
    c_flank_pwm_data, orient="index", columns=["Pos1", "Pos2", "Pos3"]
)

n_flank_pwm_data = {
    "A": [0.672938, 0.494511, 0.290216],
    "C": [-5.464582, -5.140732, -5.112201],
    "D": [0.856850, 0.732683, 0.398964],
    "E": [0.692225, 0.660452, 0.746641],
    "F": [-1.024461, -1.529751, -0.687119],
    "G": [0.873872, 0.746332, 0.630604],
    "H": [-1.386627, -1.212169, -1.011943],
    "I": [0.138351, -0.461093, 0.095419],
    "K": [0.801095, 0.639492, 0.847272],
    "L": [0.562420, -0.162599, 0.201511],
    "M": [-2.230132, -2.754557, -2.397489],
    "N": [-0.198452, -0.099572, -0.214853],
    "P": [-1.491966, 1.405721, 0.532710],
    "Q": [-0.622442, -0.151155, -0.006518],
    "R": [0.216375, 0.217991, 0.545768],
    "S": [0.623057, 0.416164, 0.236042],
    "T": [0.160517, 0.011151, -0.111494],
    "V": [0.523780, 0.239183, 0.330202],
    "W": [-3.276340, -4.050898, -2.959115],
    "Y": [-1.060520, -1.689795, -0.674071],
    # Add a pseudo-entry for 'X' with 0 contribution
    "X": [0.0, 0.0, 0.0],
}

n_flank_pwm = pd.DataFrame.from_dict(
    n_flank_pwm_data, orient="index", columns=["Pos1", "Pos2", "Pos3"]
)


class PWMFeatureGenerator(BaseFeatureGenerator):
    """
    Generates PWM (Position Weight Matrix) features for peptides based on specified MHC alleles.

    This generator calculates PWM scores for each peptide against the provided MHC class I or II allele PWMs.
    """

    CLASS_II_CORE_LENGTH = 9
    DEFAULT_PWM_PATH = os.path.join(os.path.dirname(__file__), "..", "PWMs")

    def __init__(
        self,
        peptides: List[str],
        alleles: List[str],
        anchors: int = 2,
        mhc_class: str = "I",
        pwm_path: Optional[Union[str, os.PathLike]] = None,
        remove_pre_nxt_aa: bool = False,
        remove_modification: bool = True,
        *args,
        **kwargs,
    ):
        """
        Initializes the PWMFeatureGenerator.

        Parameters:
            peptides (List[str]): Series of peptide sequences.
            alleles (List[str]): List of MHC allele names (e.g., ['HLA-A01:01', 'HLA-B07:02']).
            mhc_class (str): MHC class, either 'I' or 'II'. Default is 'I'.
            pwm_path (Optional[Union[str, os.PathLike]]): Custom path to PWM files. Defaults to '../../data/PWMs'.
            remove_pre_nxt_aa (bool): Whether to include the previous and next amino acids in peptides.
                If True, remove them. Default is True.
            remove_modification (bool): Whether to include modifications in peptides.
                If True, remove them. Default is True.
        """
        self.peptides = pd.Series(peptides)
        self.alleles = alleles
        self.mhc_class = mhc_class.upper()
        if self.mhc_class not in {"I", "II"}:
            raise ValueError("MHC class must be 'I' or 'II'.")
        self.pwm_path = pwm_path if pwm_path else PWMFeatureGenerator.DEFAULT_PWM_PATH
        logger.info(f"PWM path: {self.pwm_path}")
        self.pwms: Dict[str, Dict[int, pd.DataFrame]] = (
            self._load_pwms()
        )  # Dict[allele, Dict[mer, pd.DataFrame]]

        # logger all the pwms
        for allele, pwms in self.pwms.items():
            for mer, pwm in pwms.items():
                logger.debug(
                    f"Loaded PWM for allele {allele}, length {mer}: {pwm.shape[1]} positions"
                )
                logger.debug(pwm)
        self.anchors = anchors
        if mhc_class == "I" and anchors > 0:
            logger.info("Number of anchors: {}".format(self.anchors))
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        logger.info(
            "Remove pre and next amino acids: {}".format(self.remove_pre_nxt_aa)
        )
        self.remove_modification = remove_modification
        logger.info("Remove modifications: {}".format(self.remove_modification))
        logger.info(
            f"Initialized PWMFeatureGenerator with {len(peptides)} peptides and alleles: {alleles}"
        )

    @property
    def id_column(self) -> List[str]:
        """
        Returns a list of input columns required for the feature generator.
        """
        return ["Peptide"]

    def _extract_trailing_numbers(self, text: str) -> str:
        """
        Extracts the trailing numbers from a string.

        Examples:
            - 'ABC123' -> '123'
            - 'ABC' -> None
            - 'L123' -> '123'
        """
        import re

        match = re.search(r"\d+$", text)
        return match.group() if match else None

    def _default_allele_pwm_files(self) -> Dict[str, Dict[int, str]]:
        """
        Constructs default PWM file paths for each allele based on MHC class.
        """
        class_path = os.path.join(self.pwm_path, self.mhc_class)
        pwm_files = {}
        for allele in self.alleles:
            pwm_files[allele] = {}
            if self.mhc_class == "I":
                allele_dir = allele.replace("*", "").replace(":", "_")
            elif self.mhc_class == "II":
                allele_dir = allele.replace("*", "").replace(":", "")
                if len(allele_dir) == 8:
                    # DRB10101 -> DRB1_0101
                    allele_dir = f"{allele_dir[:4]}_{allele_dir[4:]}"
            allele_dir_path = os.path.join(class_path, allele_dir)
            pwm_files_list = utils.list_all_files_in_directory(allele_dir_path)
            logger.info(
                f"Searching for PWM files for allele {allele} in {allele_dir_path}"
            )
            logger.info(f"Found PWM files for allele {allele}: {pwm_files_list}")
            if self.mhc_class == "I":
                for pwm_file in pwm_files_list:
                    # assume the trailing numbers in the file name indicate the mer length
                    mer = int(self._extract_trailing_numbers(str(pwm_file)))
                    if mer is not None:
                        pwm_files[allele][mer] = pwm_file
                    else:
                        logger.error(
                            f"Mer length not found in PWM file name: {pwm_file}. Assuming the tailing number is the mer length."
                        )
                        raise ValueError(
                            f"Mer length not found in PWM file name: {pwm_file}"
                        )
            elif self.mhc_class == "II":
                # class II PWMs are fixed length: 9, which indicates the core length
                if len(pwm_files_list) == 1:
                    pwm_files[allele][PWMFeatureGenerator.CLASS_II_CORE_LENGTH] = (
                        pwm_files_list[0]
                    )
                else:
                    logger.error(
                        f"Expected 1 PWM file for allele {allele}, found {len(pwm_files_list)}: {pwm_files_list}"
                    )
                    raise ValueError(
                        f"Expected 1 PWM file for allele {allele}, found {len(pwm_files_list)}: {pwm_files_list}"
                    )

        logger.info(f"Default PWM file paths set for alleles: {self.alleles}")
        return pwm_files

    def _most_conserved_postions(self, pwm: pd.DataFrame, n: int = 2) -> List[int]:
        """
        Returns the n most conserved positions in the PWM.
        In our study, we only use the anchor score for class I MHC.
        """
        if n > pwm.shape[1]:
            raise ValueError(
                f"Number of positions to return ({n}) exceeds the PWM length ({pwm.shape[1]})."
            )
        pfm = pwm.apply(lambda x: 2**x)
        entropy = -1 * (pfm * np.log2(pfm)).sum(axis=0)
        return entropy.nsmallest(n).index.tolist()

    def _load_pwms(self) -> Dict[str, Dict[int, pd.DataFrame]]:
        """
        Loads PWMs for each allele from the constructed file paths.
        """
        pwms = {}
        allele_pwm_files = self._default_allele_pwm_files()
        for allele, mer_files in allele_pwm_files.items():
            pwms[allele] = {}
            for mer, file_path in mer_files.items():
                if not os.path.exists(file_path):
                    logger.error(f"PWM file not found: {file_path}")
                    raise FileNotFoundError(f"PWM file not found: {file_path}")
                try:
                    pwm = pd.read_csv(
                        file_path, delim_whitespace=True, header=None, index_col=0
                    )
                    pwm.columns = [f"{pos+1}" for pos in range(pwm.shape[1])]
                    pwms[allele][mer] = pwm
                    logger.info(
                        f"Loaded PWM for allele {allele}, length {mer} from {file_path}"
                    )
                except Exception as e:
                    logger.error(f"Error loading PWM file {file_path}: {e}")
                    raise e
        return pwms

    def _cal_PWM_score_I(self, peptide: str, allele: str) -> Optional[float]:
        """
        Calculates PWM scores for MHC class I.
        If length is out of range, return pd.NA.
        """
        peptide_len = len(peptide)
        min_mer = min(self.pwms[allele].keys())
        max_mer = max(self.pwms[allele].keys())

        if peptide_len < min_mer or peptide_len > max_mer:
            return pd.NA
        else:
            pwm = self.pwms[allele][peptide_len]
            try:
                score = sum(pwm.loc[aa, str(i + 1)] for i, aa in enumerate(peptide))
            except KeyError as e:
                logger.error(
                    f"Residue '{e}' not found in PWM for allele={allele}, "
                    f"peptide={peptide}, length={peptide_len}"
                )
                raise ValueError(f"Residue '{e}' not found in PWM.")
            return score

    def _cal_PWM_score_II(
        self, peptide: str, allele: str
    ) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """
        Calculates PWM scores for MHC class II using a sliding 9-mer window.

        We:
          1) Slide over all possible 9-mer windows to find the highest core PWM score.
          2) Once the best core is found, extract up to 3 AA on each flank (N-flank and C-flank).
          3) If the flank has fewer than 3 residues, pad with 'X'.
          4) Score each flank with n_flank_pwm and c_flank_pwm and return all three scores.

        Returns:
            (core_score, n_flank_score, c_flank_score)
            or (NaN, NaN, NaN) if peptide has length < 9.
        """
        core_len = PWMFeatureGenerator.CLASS_II_CORE_LENGTH
        pwm = self.pwms[allele][core_len]

        if len(peptide) < core_len:
            # Return NaNs for the scores if peptide too short
            return (pd.NA, pd.NA, pd.NA)

        best_score = None
        best_core = None
        best_core_start_idx = 0

        # Slide over possible 9-mer windows
        for start_idx in range(len(peptide) - core_len + 1):
            subpeptide_9 = peptide[start_idx : start_idx + core_len]
            try:
                tmp_score = sum(
                    float(pwm.loc[aa, str(i + 1)]) for i, aa in enumerate(subpeptide_9)
                )
            except KeyError as e:
                logger.error(
                    f"Residue '{e}' not found in PWM for allele={allele}, "
                    f"peptide={peptide} (subpep={subpeptide_9})."
                )
                raise ValueError(f"Residue '{e}' not found in PWM.")

            if (best_score is None) or (tmp_score > best_score):
                best_score = tmp_score
                best_core = subpeptide_9
                best_core_start_idx = start_idx

        logger.debug(
            f"Peptide: {peptide}, best core: {best_core}, "
            f"start_idx: {best_core_start_idx}, score: {best_score}"
        )

        # Compute N-flank and C-flank (up to 3 residues)
        n_flank_start = max(0, best_core_start_idx - 3)
        n_flank_end = best_core_start_idx
        c_flank_start = best_core_start_idx + core_len
        c_flank_end = min(len(peptide), c_flank_start + 3)

        n_flank_seq = peptide[n_flank_start:n_flank_end]
        c_flank_seq = peptide[c_flank_start:c_flank_end]

        # Pad with 'X' if flank < 3 residues
        n_flank_seq = (("X" * (3 - len(n_flank_seq))) + n_flank_seq)[
            -3:
        ]  # ensure length = 3
        c_flank_seq = (c_flank_seq + ("X" * (3 - len(c_flank_seq))))[
            :3
        ]  # ensure length = 3

        n_flank_score = 0.0
        c_flank_score = 0.0
        for i, aa in enumerate(n_flank_seq):
            try:
                n_flank_score += n_flank_pwm.loc[aa, f"Pos{i+1}"]
            except KeyError:
                # 'X' or unknown -> 0
                pass

        for i, aa in enumerate(c_flank_seq):
            try:
                c_flank_score += c_flank_pwm.loc[aa, f"Pos{i+1}"]
            except KeyError:
                # 'X' or unknown -> 0
                pass

        logger.debug(
            f"n_flank: {n_flank_seq}, c_flank: {c_flank_seq}, "
            f"n_flank_score: {n_flank_score}, c_flank_score: {c_flank_score}"
        )

        return (best_score, n_flank_score, c_flank_score)

    def _cal_PWM_score(
        self, peptide: str, allele: str
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculates PWM scores for a single peptide across all applicable mer lengths for a given allele.

        For MHC class I:
            Returns a single float (or pd.NA).
        For MHC class II:
            Returns a tuple of (core_score, n_flank_score, c_flank_score) or (pd.NA, pd.NA, pd.NA).
        """
        if self.mhc_class == "I":
            return self._cal_PWM_score_I(peptide, allele)
        else:
            return self._cal_PWM_score_II(peptide, allele)

    def _cal_anchor_score(
        self, peptide: str, allele: str, anchor_dict: Dict[int, List[int]]
    ) -> Optional[float]:
        """
        Calculate anchor score for a single peptide across all applicable mer lengths for a given allele and given positions.

        Parameters:
            peptide (str): The peptide sequence.
            allele (str): The MHC allele to score against.
            anchor_dict (Dict[int, List[int]]): Dictionary containing the most conserved positions for each mer length.

        Returns:
            Optional[float]: Anchor score for the peptide against the allele's PWM, or NaN if out of range.
        """
        peptide_len = len(peptide)
        if self.mhc_class == "I":
            min_mer = min(self.pwms[allele].keys())
            max_mer = max(self.pwms[allele].keys())
            if peptide_len < min_mer or peptide_len > max_mer:
                return pd.NA
            else:
                pwm = self.pwms[allele].get(peptide_len)
                score = 0
                anchors = anchor_dict[peptide_len]
                try:
                    for anchor in anchors:
                        anchor = int(anchor)
                        score += pwm.loc[peptide[anchor - 1], str(anchor)]
                except KeyError as e:
                    logger.error(
                        f"Position '{e}' not found in PWM for allele {allele} and peptide {peptide} with length {peptide_len}."
                    )
                    raise ValueError(f"Position '{e}' not found in PWM.")
                return score
        elif self.mhc_class == "II":
            logger.warning("Anchor score calculation not implemented for MHC class II.")
            return None

    def set_pwms(self, pwms: Dict[str, Dict[int, pd.DataFrame]]):
        """
        Sets PWMs directly, allowing for custom PWMs to be provided.
        """
        self.pwms = pwms
        logger.info(f"Set custom PWMs for alleles: {list(pwms.keys())}")

    def generate_features(self) -> pd.DataFrame:
        """
        Generates PWM features for all peptides across specified alleles.

        For MHC class I:
            - 'PWM_Score_{allele}' and optionally 'Anchor_Score_{allele}' columns.
        For MHC class II:
            - 'PWM_Score_{allele}' (core 9-mer),
            - 'N_Flank_PWM_Score_{allele}',
            - 'C_Flank_PWM_Score_{allele}' columns.
        """
        features_df = pd.DataFrame(self.peptides, columns=["Peptide"])
        features_df["clean_peptide"] = features_df["Peptide"]
        if self.remove_pre_nxt_aa:
            features_df["clean_peptide"] = features_df["Peptide"].apply(
                utils.remove_pre_and_nxt_aa
            )
        if self.remove_modification:
            features_df["clean_peptide"] = features_df["clean_peptide"].apply(
                utils.remove_modifications
            )

        # Convert nonstandard amino acids: U -> C
        features_df["clean_peptide"] = features_df["clean_peptide"].apply(
            lambda x: x.replace("U", "C")
        )

        for allele in self.alleles:
            logger.info(
                f"Generating PWM scores for allele: {allele}, total peptides: {len(features_df)}"
            )

            if self.mhc_class == "I":
                # Class I returns a single score
                features_df[f"PWM_Score_{allele}"] = features_df["clean_peptide"].apply(
                    lambda peptide: self._cal_PWM_score(peptide, allele)
                )
                na_count = features_df[f"PWM_Score_{allele}"].isna().sum()
                logger.info(
                    f"Missing PWM scores for {na_count} peptides. Using median for imputation."
                )
                features_df.fillna(
                    {
                        f"PWM_Score_{allele}": features_df[
                            f"PWM_Score_{allele}"
                        ].median()
                    },
                    inplace=True,
                )

                if self.anchors != 0:
                    logger.info(
                        f"Generating anchor scores for allele: {allele}, total peptides: {len(features_df)}"
                    )
                    anchor_dict = {}
                    min_mer = min(self.pwms[allele].keys())
                    max_mer = max(self.pwms[allele].keys())
                    for mer_len in range(min_mer, max_mer + 1):
                        anchor_dict[mer_len] = self._most_conserved_postions(
                            self.pwms[allele][mer_len], self.anchors
                        )
                    logger.info(
                        f"Most conserved positions for allele {allele}: {anchor_dict}"
                    )
                    features_df[f"Anchor_Score_{allele}"] = features_df[
                        "clean_peptide"
                    ].apply(
                        lambda peptide: self._cal_anchor_score(
                            peptide, allele, anchor_dict
                        )
                    )
                    na_count = features_df[f"Anchor_Score_{allele}"].isna().sum()
                    logger.info(
                        f"Missing anchor scores for {na_count} peptides. Using median for imputation."
                    )
                    features_df.fillna(
                        {
                            f"Anchor_Score_{allele}": features_df[
                                f"Anchor_Score_{allele}"
                            ].median()
                        },
                        inplace=True,
                    )

            else:
                # Class II returns (core_score, n_flank_score, c_flank_score)
                features_df[
                    [
                        f"PWM_Score_{allele}",
                        f"N_Flank_PWM_Score_{allele}",
                        f"C_Flank_PWM_Score_{allele}",
                    ]
                ] = features_df["clean_peptide"].apply(
                    lambda pep: pd.Series(self._cal_PWM_score(pep, allele))
                )

                # Impute missing (NaN) with medians for each new column
                for col in [
                    f"PWM_Score_{allele}",
                    f"N_Flank_PWM_Score_{allele}",
                    f"C_Flank_PWM_Score_{allele}",
                ]:
                    na_count = features_df[col].isna().sum()
                    logger.info(
                        f"Missing {col} for {na_count} peptides. Using median for imputation."
                    )
                    median_val = features_df[col].median()
                    features_df[col].fillna(median_val, inplace=True)

        features_df.drop(columns=["clean_peptide"], inplace=True)

        return features_df

    @property
    def feature_columns(self) -> List[str]:
        """
        Returns a list of feature names generated by the feature generator.
        """
        if self.mhc_class == "I":
            feature_columns = [f"PWM_Score_{allele}" for allele in self.alleles]
            if self.anchors != 0:
                anchor_columns = [f"Anchor_Score_{allele}" for allele in self.alleles]
                feature_columns.extend(anchor_columns)
            return feature_columns
        else:
            # Class II
            feature_columns = []
            for allele in self.alleles:
                feature_columns.append(f"PWM_Score_{allele}")
                feature_columns.append(f"N_Flank_PWM_Score_{allele}")
                feature_columns.append(f"C_Flank_PWM_Score_{allele}")
            return feature_columns
