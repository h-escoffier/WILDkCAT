import os
import time
import socket
import logging
import pandas as pd
from Bio import Align, Entrez
from dotenv import load_dotenv
from functools import lru_cache
from urllib.error import HTTPError, URLError

from ..api.uniprot_api import convert_uniprot_to_sequence   


load_dotenv()


def closest_enz(kcat_dict, api_output) -> pd.DataFrame:
    """
    Retrieve and ranks the enzymes sequences closest to the sequence of the target enzyme based on the percentage of identity.
    If the reference UniProt ID is missing, invalid, or the sequence cannot be retrieved, the function returns the input DataFrame with "id_perc" set to None.
    
    Parameters:    
        kcat_dict (dict): Dictionary containing at least the key 'uniprot_model' with the reference UniProt ID.
        api_output (pd.DataFrame): DataFrame containing a column "UniProtKB_AC" with UniProt IDs to compare against.
    
    Returns:
        pd.DataFrame: A copy of `api_output` with an added "id_perc" column (identity percentage). 
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

    ref_uniprot_id = kcat_dict.get('catalytic_enzyme')
    if pd.isna(ref_uniprot_id) or (";" in str(ref_uniprot_id)):
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


Entrez.email = os.getenv("ENTREZ_EMAIL")


@lru_cache(maxsize=None)
def _fetch_taxonomy(species_name):
    """Fetch and cache the taxonomic lineage of a species from NCBI."""
    for attempt in range(3):
        try:
            with Entrez.esearch(db="taxonomy", term=species_name) as handle:
                record = Entrez.read(handle)
            if not record["IdList"]:
                return ()
            tax_id = record["IdList"][0]
            with Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml") as handle:
                records = Entrez.read(handle, validate=False)
            if not records:
                return ()
            lineage = tuple(
                taxon["ScientificName"] for taxon in records[0]["LineageEx"]
            ) + (records[0]["ScientificName"],)
            return lineage
        except (HTTPError, URLError, socket.error, RuntimeError) as e: # TODO: Put in log file ? 
            print(f"[Retry {attempt+1}/3] Taxonomy fetch failed for '{species_name}': {e}")
            time.sleep(3)
        except Exception as e:
            print(f"[Fatal] '{species_name}': {type(e).__name__}: {e}")
            return ()
    return ()


@lru_cache(maxsize=None)
def _calculate_taxonomy_score(ref_organism, target_organism):
    """
    Calculate a taxonomy distance score between reference and target organisms.
        
    Parameters: 
        ref_organism (str): The reference organism's name.
        target_organism (str): The target organism's name.

    Returns:
        int: distance between reference and target organisms (0 = identical species, higher = more distant).
    """
    ref_lineage = _fetch_taxonomy(ref_organism)
    target_lineage = _fetch_taxonomy(target_organism)

    if not target_lineage:
        return len(ref_lineage) + 1  # Penalize missing taxonomy

    ref_set = set(ref_lineage)
    similarity = 0
    for taxon in target_lineage:
        if taxon in ref_set:
            similarity += 1
        else:
            break
    return len(ref_lineage) - similarity


def closest_taxonomy(general_criteria, api_output) -> pd.DataFrame:
    """
    Retrieve and ranks the organisms based on their taxonomic similarity to the reference organism.
    
    Parameters:    
        general_criteria (dict): Dictionary containing at least the key 'organism' with the reference organism.
        api_output (pd.DataFrame): DataFrame containing a column "Organism". 
    
    Returns:
        pd.DataFrame: A copy of `api_output` with an added "organism_score" column.
    """
    ref_organism = general_criteria["Organism"]
    api_output = api_output.copy()

    unique_organisms = api_output["Organism"].unique()
    score_map = {
        organism: _calculate_taxonomy_score(ref_organism, organism)
        for organism in unique_organisms
    }

    api_output["organism_score"] = api_output["Organism"].map(score_map)
    return api_output