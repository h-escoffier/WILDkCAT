import requests
import logging
import pandas as pd 
from tqdm import tqdm
import re
from io import StringIO
from functools import lru_cache

from kcatmatchmod.api.api_utilities import retry_api, safe_requests_get


# TODO: Integrate safe requests and retry_api decorators in the functions below


# --- API --- 


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


@lru_cache(maxsize=None)
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
        logging.warning('%s: Failed to retrieve SID for KEGG compound ID %s' % (kegg_compound_id))
        return None
    cid = convert_sid_to_cid(sid)
    if cid is None:
        logging.warning('%s: Failed to retrieve CID for KEGG compound ID %s' % (kegg_compound_id))
        return None
    smiles = convert_cid_to_smiles(cid)
    if smiles is None:
        logging.warning('%s: Failed to retrieve SMILES for KEGG compound ID %s' % (kegg_compound_id))
        return None
    return smiles
    

# --- Retrieve Sequences from UniProtID --- 

@lru_cache(maxsize=None)
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


# --- Create CataPro input file ---


def create_catapro_input_file(kcat_df, output_path):
    """
    TODO: 
    """
    catapro_input = []

    for _, row in tqdm(kcat_df.iterrows(), total=len(kcat_df), desc="Generating CataPro input"):
        uniprot = row['uniprot_model']
        ec_code = row['ec_code']

        # If multiple UniProt IDs continue TODO: Find a way to handle this 
        if len(uniprot.split(';')) > 1:        
            logging.warning(f"Multiple UniProt IDs found for {ec_code}: {uniprot}.")
            continue

        sequence = convert_uniprot_to_sequence(uniprot) 
        if sequence is None:
            continue
        
        smiles_list = []
        
        for kegg_compound_id in row['substrates_kegg'].split(';'):
            smiles = convert_kegg_to_smiles(kegg_compound_id)
            if smiles is not None:
                smiles_list.append(smiles[0]) # If multiple SMILES, take the first one TODO: Handle this case
        
        if len(smiles_list) > 0:
            for smiles in smiles_list:
                catapro_input.append({
                    "Enzyme_id": uniprot,
                    "type": "wild",  # TODO: I think it should be always 'wild' ? 
                    "sequence": sequence,
                    "smiles": smiles
                })

    # Generate CataPro input file
    catapro_input_df = pd.DataFrame(catapro_input)
    catapro_input_df.to_csv(output_path, sep=',', index=True)

    return catapro_input_df


# --- Integrate CataPro predictions into kcat file ---

def integrate_catapro_predictions(kcat_file_path, catapro_predictions_path, output_path):
    pass 


# --- Main ---


def run_catapro(kcat_file_path, limit_matching_score, output_path, report=True):
    """
    TODO 
    """
    # Read the kcat file
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')

    # Subset rows with no values or matching score above the limit
    # kcat_df = kcat_df[(kcat_df['matching_score'] < limit_matching_score) | (kcat_df['matching_score'].isnull())]
    # Drop rows with no UniProt ID or no substrates_kegg
    kcat_df = kcat_df[kcat_df['uniprot_model'].notnull() & kcat_df['substrates_kegg'].notnull()]
    
    # Generate CataPro input file
    catapro_input = create_catapro_input_file(kcat_df, output_path)

    # Run CataPro predictions
    # Integrate CataPro predictions into kcat file
    # Save the output file
    # Generate report if required


if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))

    # Test : Main function
    run_catapro("output/ecoli_kcat_sabio.tsv", 9, "in_progress/catapro_input.csv")