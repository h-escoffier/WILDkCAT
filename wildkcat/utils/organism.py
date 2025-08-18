import pandas as pd
import numpy as np
from Bio import Align


from wildkcat.api.uniprot_api import convert_uniprot_to_sequence   


def closest_org(kcat_dict, api_output):
    """
    Retrieve and ranks the enzymes sequences closest to the sequence of the target enzyme based on the percentage of identity.
    If the reference UniProt ID is missing, invalid, or the sequence cannot be retrieved, the function returns the input DataFrame with "id_perc" set to None.
    
    Parameters:    
        kcat_dict (dict): Dictionary containing at least the key 'uniprot_model' with the reference UniProt ID.
        api_output (pd.DataFrame): DataFrame containing a column "UniProtKB_AC" with UniProt IDs to compare against.
    
    Returns:
        pd.DataFrame: A copy of `api_output` with an added "id_perc" column (identity percentage), sorted by identity in descending order.
    """

    def _calculate_identity(seq_ref, seq_db):
        """
        Returns the percentage of identical characters between two sequences.
        Adapted from https://gist.github.com/JoaoRodrigues/8c2f7d2fc5ae38fc9cb2 

        Parameters: 
            seq_ref (str): The reference sequence.
            seq_db (str): The sequence to compare against.

        Returns: 
            float: The percentage of identical characters between the two sequences.
        """
        matches = [a == b for a, b in zip(seq_ref, seq_db)]
        return (100 * sum(matches)) / len(seq_ref)

    ref_uniprot_id = kcat_dict.get('uniprot_model')
    if pd.isna(ref_uniprot_id) or ';' in str(ref_uniprot_id):
        api_output = api_output.copy()
        api_output["id_perc"] = None
        return api_output
    
    ref_seq = convert_uniprot_to_sequence(ref_uniprot_id)
    if ref_seq is None:
        api_output = api_output.copy()
        api_output["id_perc"] = None
        return api_output

    aligner = Align.PairwiseAligner()
    identity_scores = []
    
    for uniprot_id in api_output["UniProtKB_AC"]:
        if pd.isna(uniprot_id):
            identity_scores.append(None)
            continue
        seq = convert_uniprot_to_sequence(uniprot_id)
        if seq is None:
            identity_scores.append(None)
            continue
        elif len(seq) == 0:
            identity_scores.append(0)
            continue

        alignments = aligner.align(ref_seq, seq)
        aligned_ref, aligned_db = alignments[0]
        id_score = _calculate_identity(aligned_ref, aligned_db)
        identity_scores.append(id_score)

    api_output = api_output.copy()
    api_output["id_perc"] = identity_scores

    return api_output