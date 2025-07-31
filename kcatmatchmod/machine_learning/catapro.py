import requests
import logging
import pandas as pd 
from tqdm import tqdm
import re
from io import StringIO
from functools import lru_cache


# ---------------------------------------------
# Retrieve SMILES from KEGG ID
# ---------------------------------------------


def convert_kegg_compound_to_sid(kegg_compound_id):
    """
    Convert the KEGG compound ID to the PubChem Substance ID (SID).

    Parameters:
        kegg_compound_id (str): KEGG compound ID.

    Returns:
        str: The PubChem SID if found, otherwise None.
    """
    url = f"https://rest.kegg.jp/conv/pubchem/compound:{kegg_compound_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None

    match = re.search(r'pubchem:\s*(\d+)', response.text)
    sid = match.group(1) if match else None
    return sid


def convert_sid_to_cid(sid):
    """
    Converts a PubChem Substance ID (SID) to the corresponding Compound ID (CID).

    Parameters:
        sid (str): PubChem Substance ID.

    Returns:
        int or None: The corresponding PubChem Compound ID (CID), or None if not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/JSON"
    response = requests.get(url)
    cid = None
    if response.status_code == 200:
        try:
            cid = response.json()['InformationList']['Information'][0]['CID'][0]
        except (KeyError, IndexError):
            cid = None
    return cid


def convert_cid_to_smiles(cid):    
    """
    Converts a PubChem Compound ID (CID) to its corresponding SMILES representation.

    Parameters:
        cid (str): PubChem Compound ID.

    Returns:
       list 
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/smiles/txt"
    try:
        response = requests.get(url)
        response.raise_for_status()
        smiles = response.text.strip().split('\n')
        return smiles
    except:
        return None


def convert_kegg_to_smiles(kegg_compound_id):
    """
    Convert the KEGG compound ID to the PubChem Compound ID (CID).

    Parameters:
        kegg_compound_id (str): KEGG compound ID.

    Returns:
        list
    """
    sid = convert_kegg_compound_to_sid(kegg_compound_id)
    if sid is None:
        return None
    cid = convert_sid_to_cid(sid)
    if cid is None:
        return None
    return convert_cid_to_smiles(cid)
    

@lru_cache(maxsize=None)
def cached_convert_kegg_to_smiles(kegg_compound_id):
    return convert_kegg_to_smiles(kegg_compound_id)


# ---------------------------------------------
# Retrieve Sequences from UniProtID 
# ---------------------------------------------

def convert_uniprot_to_sequence(uniprot_id):
    """
    Convert a UniProt accession ID to its corresponding amino acid sequence.

    Parameters:
        uniprot_id (str): The UniProt accession ID.

    Returns:
        str: The amino acid sequence, or None if not found.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        fasta = response.text
        lines = fasta.splitlines()
        sequence = ''.join(lines[1:])  # Skip the header
        return sequence
    else:
        logging.warning('%s: Failed to retrieve sequence for UniProt ID %s' % (uniprot_id))
        return None


@lru_cache(maxsize=None)
def cached_convert_uniprot_to_sequence(uniprot_id):
    return convert_uniprot_to_sequence(uniprot_id)


# ---------------------------------------------
# Create CataPro input file 
# ---------------------------------------------


def create_catapro_input_file(kcat_file_path, output_path):
    pass 


# ---------------------------------------------
# Integrate CataPro predictions into kcat file 
# ---------------------------------------------

def integrate_catapro_predictions(kcat_file_path, catapro_predictions_path, output_path):
    pass 


# ---------------------------------------------
# Main functions
# ---------------------------------------------

def run_catapro(kcat_file_path, output_path, report=True):
    """
    TODO 
    """
    # Load kcat file 
    # Identify the kcat values to predict
    # Retrieve SMILES and sequences
    # Create CataPro input file
    # Run CataPro predictions
    # Integrate CataPro predictions into kcat file
    # Save the output file
    # Generate report if required
    pass


if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    print(convert_uniprot_to_sequence("P0A796"))

    # Test : Main function
    # run_catapro("input/kcat_file.tsv", "output/catapro_input.tsv