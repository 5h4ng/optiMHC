import logging
from typing import List, Optional, Tuple, Union
import pandas as pd
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


def read_pin(
    pin_files: Union[str, List[str], pd.DataFrame],
    score_column: Optional[str] = None,
    retention_time_column: Optional[str] = None,
    *args,
    **kwargs
) -> PsmContainer: 
    """
    Reads a Percolator INput (PIN) file into a PsmContainer object.

    Parameters:
        pin_file (Union[str, List[str]], pd.DataFrame): The file path to the PIN file or a list of file paths.
        score_column (Optional[str]): The column containing the Percolator score.
        retention_time_column (Optional[str]): The column containing the retention time.
    
    Returns:
        PsmContainer: A PsmContainer object containing the PSM data.
    """
    logger.info("Reading PIN file(s) into PsmContainer.")
    if isinstance(pin_files, str):
        pin_files = [pin_files]

    pin_df = pd.concat(
        [_read_single_pin_as_df(pin_file) for pin_file in pin_files]
    )
    logger.info(f"Read {len(pin_df)} PSMs from {len(pin_files)} PIN files.")
    logger.debug(pin_df.head())

    def find_required_columns(col: str, columns: List[str]):
        if col not in columns:
            raise ValueError(f"Column '{col}' not found in PSM data.")
        return col
    
    # non-feature columns
    label = find_required_columns('Label', pin_df.columns)
    scan = find_required_columns('ScanNr', pin_df.columns)
    specid = find_required_columns('SpecId', pin_df.columns)
    peptide = find_required_columns('Peptide', pin_df.columns)
    protein = find_required_columns('Proteins', pin_df.columns)

    # Comet: P2PI20160713_pilling_C1RA2_BB72_P1_31_3_1
    # - unique_id: P2PI20160713_pilling_C1RA2_BB72_P1_31_3
    # - hit_rank: 1

    # MSFragger: P2PI20160713_pilling_C1RA2_BB72_P1.3104.3104.2_1
    # - unique_id: P2PI20160713_pilling_C1RA2_BB72_P1.3104.3104.2
    # - hit_rank: 1

    def parse_specid(specid: str) -> Tuple[str, int]:
        parts = specid.rsplit('_', 1)
        unique_id = parts[0]
        hit_rank = int(parts[1])
        return unique_id, hit_rank
    
    hit_rank = 'rank'
    if 'rank' in pin_df.columns: # MSFragger
        pin_df[specid], _ = zip(*pin_df[specid].apply(parse_specid))
    else:
        pin_df[specid], pin_df['rank'] = zip(*pin_df[specid].apply(parse_specid))

    # feature columns: columns that are not non-feature columns
    retention_time_column = find_required_columns(retention_time_column, pin_df.columns) if retention_time_column else None
    feature_columns = [col for col in pin_df.columns if col not in [label, scan, specid, peptide, protein]]

    logger.info(f"Columns: label={label}, scan={scan}, specid={specid}, peptide={peptide}, protein={protein}, hit_rank={hit_rank}, score={score_column}, retention_time={retention_time_column}, features={feature_columns}")

    # convert the type
    pin_df[scan] = pin_df[scan].astype(str)
    pin_df[specid] = pin_df[specid].astype(str)
    pin_df[peptide] = pin_df[peptide].astype(str)
    pin_df[protein] = pin_df[protein].astype(str)
    pin_df[hit_rank] = pin_df[hit_rank].astype(int)
    if retention_time_column:
        pin_df[retention_time_column] = pin_df[retention_time_column].astype(float)
    for col in feature_columns:
        pin_df[col] = pin_df[col].astype(float)
    
    # label = 1 for target, -1 for decoy. Convert to Boolean.
    pin_df[label] = (pin_df[label] == '1')
    rescoring_features = {'Original': feature_columns}

    return PsmContainer(
        psms=pin_df,
        label_column=label,
        scan_column=scan,
        spectra_column=specid,
        ms_data_file_column=None,
        peptide_column=peptide,
        protein_column=protein,
        rescoring_features=rescoring_features,
        hit_rank_column=hit_rank,
        retention_time_column=retention_time_column
    )
  
  
def _read_single_pin_as_df(pin_file: str) -> pd.DataFrame:
    """
    Proteins column in PIN file is a tab-separated list of proteins.
    This function reads the PIN file and store the proteins in one column of dataframe.

    Parameters:
        pin_file (str): The file path to the PIN file.
    
    Returns:
        pd.DataFrame: A DataFrame containing the PSM data.
    """
    logger.info(f"Reading PIN file: {pin_file}")
    with open(pin_file, 'r') as f:
        header = f.readline().strip().split("\t")
        header_len = len(header)
        data = []
        for line in f:
            parts = line.strip().split("\t")
            proteins_column_num  = len(parts) - header_len + 1
            proteins = '\t'.join(parts[-proteins_column_num:])
            data.append(parts[:len(parts) - proteins_column_num] + [proteins])
    df = pd.DataFrame(data, columns=header)
    return df

