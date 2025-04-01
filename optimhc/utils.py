# utils.py

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Union
from logging import getLogger

logger = getLogger(__name__)


def convert_pfm_to_pwm(pfm_filename, pseudocount=0.8, background_freqs=None):
    """
    Converts a Position Frequency Matrix (PFM) file to a Position Weight Matrix (PWM).

    Parameters:
        pfm_filename (str): The file path to the PFM file.
        pseudocount (float): The pseudocount to add to the PFM to avoid zero probabilities.
        background_freqs (dict): A dictionary containing the background frequencies for each nucleotide.

    Returns:
        pd.DataFrame: A DataFrame representation of the PWM.
    """
    # Default background frequencies if not provided
    if background_freqs is None:
        background_freqs = 1 / 20

    pfm = pd.read_csv(pfm_filename, sep="\t", header=None, index_col=0)
    pfm.drop(pfm.columns[-1], axis=1, inplace=True)  # Drop any extraneous columns
    pfm_pseudo = pfm + pseudocount
    ppm = pfm_pseudo.div(pfm_pseudo.sum(axis=0), axis=1)
    pwm = np.log2(ppm.div(list(background_freqs.values())))

    return pwm


def remove_pre_and_nxt_aa(peptide: str) -> str:
    """
    Removes the pre and next amino acids from a peptide sequence.

    Example:
        - "A.AB[15.99]CDEFGHI.K" -> "ABCDEFGHI"

    Parameters:
        peptide (str): The peptide sequence.

    Returns:
        str: The peptide sequence with pre and next amino acids removed.
    """
    import re
    return re.sub(r'^[^.]*\.|\.[^.]*$', '', peptide)


def remove_modifications(peptide: str, keep_modification = None) -> str:
    """
    Removes modifications from a peptide sequence, with an option to keep specific modifications.

    Example:
        - "A.AB[15.99]CDEFGHI.K" -> "A.ABCDEFGHI.K"
        - "A.AB[UNIMOD:4]C[15.99]DEFGHI.K", keep_modification='UNIMOD:4' 
            -> "A.AB[UNIMOD:4]CDEFGHI.K"
    
    Parameters:
        peptide (str): The peptide sequence.
        keep_modification (str or list, optional): The modification(s) to keep.

    Returns:
        str: The peptide sequence with modifications removed or kept.
    """
    import re
    if keep_modification is None:
        return re.sub(r'\[.*?\]', '', peptide)
    else:
        if not isinstance(keep_modification, list):
            keep_modification = [keep_modification]
            
        def replace_mod(match):
            mod = match.group(0)
            if any(keep in mod for keep in keep_modification):
                return mod
            return ''
        return re.sub(r'\[.*?\]', replace_mod, peptide)
        

def preprocess_peptide(peptide: str) -> str:
    """
    Preprocesses the peptide sequence by removing flanking regions and modifications.

    Example:
        - "A.AB[15.99]CDEFGHI.K" -> "ABCDEFGHI"
        - "B.AC[Carbamidomethyl]DE.F" -> "ACDE"

    Parameters:
        peptide (str): Original peptide sequence.

    Returns:
        str: Cleaned peptide sequence.
    """
    peptide = remove_pre_and_nxt_aa(peptide)
    peptide = remove_modifications(peptide)
    return peptide


def list_all_files_in_directory(directory_path: str) -> List[str]:
    """
    Retrieve all files in the specified directory and return a list of file paths.

    Parameters:
        directory_path (str): The path to the directory.
        
    Returns:
        List[str]: A list of file paths.
    """
    path = Path(directory_path)
    file_list = [str(file) for file in path.rglob('*') if file.is_file()]
    return file_list


def extract_unimod_from_peptidoform(peptide: str, mod_dict: dict) -> tuple:
    """
    Converts a modified peptide sequence into DeepLC format.
    
    Parameters:
        peptide (str): The input peptide sequence with modifications in brackets.
                       Example: 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
        mod_dict (dict): A dictionary mapping modification names (in peptide) to 
                         corresponding Unimod names.
                         Example: {'15.9949': 'Oxidation', '42.0106': 'Acetyl'}
    
    Returns:
        tuple: A tuple containing:
            - seq (str): The unmodified peptide sequence.
            - modifications (str): A string of modifications formatted as 
                                   `position|UnimodName`, separated by pipes `|`.
                                   Example: '12|Oxidation|18|Acetyl'
    """
    clean_sequence = ""
    modifications = []
    current_position = 0  
    i = 0  
    while i < len(peptide):
        if peptide[i] == '[':  
            end = peptide.find(']', i)
            if end == -1:
                raise ValueError(f"Invalid modification format in {peptide}: missing closing bracket.")
            mod_name = peptide[i+1:end]
            if mod_name not in mod_dict:
                raise ValueError(f"Modification '{mod_name}' not found in the dictionary.")
            modifications.append(f"{current_position}|{mod_dict[mod_name]}")
            i = end + 1
        else:
            clean_sequence += peptide[i]
            current_position += 1
            i += 1

    modification_str = "|".join(modifications)   
    logger.debug(f'Original peptide: {peptide}')
    logger.debug(f"Output clean_sequence: {clean_sequence}")
    logger.debug(f"Output modifications: {modification_str}")
    return clean_sequence, modification_str

def convert_to_unimod_format(peptide: str, mod_dict: dict) -> str:
    """
    Converts a modified peptide sequence into Unimod format.
    
    Parameters:
        peptide (str): The input peptide sequence with modifications in brackets.
                       Example: 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
        mod_dict (dict): A dictionary mapping modification names (in peptide) to 
                         corresponding Unimod names.
                         Example: {'15.9949': 'UNIMOD:4', '42.0106': 'UNIMOD:1'}
    Returns:
        str: The peptide sequence formatted for Unimod.
             Example: 'AANDAGYFNDEM[UNIMOD:4]APIEVK[UNIMOD:1]TK'
    """
    res = peptide
    for key, value in mod_dict.items():
        res = res.replace(f'[{key}]', f'[{value}]')
    logger.debug(f'Original peptide: {peptide}')
    logger.debug(f"Output peptide: {res}")
    return res

