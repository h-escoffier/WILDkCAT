import requests
import logging
import pandas as pd 
from tqdm import tqdm
import re
from io import StringIO
from functools import lru_cache

from kcatmatchmod.api.api_utilities import retry_api, safe_requests_get
from kcatmatchmod.api.uniprot_api import get_cofactor
from kcatmatchmod.api.brenda_api import convert_uniprot_to_sequence


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
    

# --- Create CataPro input file ---


def create_catapro_input_file(kcat_df, output_path):
    """
    Generate CataPro input file and a mapping of substrate KEGG IDs to SMILES.

    Returns:
        catapro_input_df (pd.DataFrame): DataFrame for CataPro input.
        substrates_to_smiles (dict): Mapping KEGG ID <-> SMILES (both directions).
    """
    catapro_input = []
    substrates_to_smiles = {}

    for _, row in tqdm(kcat_df.iterrows(), total=len(kcat_df), desc="Generating CataPro input"):
        uniprot = row['uniprot_model']
        ec_code = row['ec_code']

        # If multiple UniProt IDs continue #TODO: Find a way to handle this 
        if len(uniprot.split(';')) > 1:        
            logging.warning(f"Multiple UniProt IDs found for {ec_code}: {uniprot}.")
            continue

        sequence = convert_uniprot_to_sequence(uniprot) 
        if sequence is None:
            continue
        
        smiles_list = []
        names = row['substrates_name'].split(';')
        kegg_ids = row['substrates_kegg'].split(';')
        
        # Get the cofactor for the EC code
        cofactor = get_cofactor(ec_code) 

        for name, kegg_compound_id in zip(names, kegg_ids):
            if name.lower() in [c.lower() for c in cofactor]:  # Skip the cofactor #TODO: Should we add a warning if no cofactor is found for a reaction?
                continue
            smiles = convert_kegg_to_smiles(kegg_compound_id)
            if smiles is not None:
                smiles_str = smiles[0]  # If multiple SMILES, take the first one #TODO: Handle this case
                smiles_list.append(smiles_str)
                substrates_to_smiles[kegg_compound_id] = smiles_str
        
        if len(smiles_list) > 0:
            for smiles in smiles_list:
                catapro_input.append({
                    "Enzyme_id": uniprot,
                    "type": "wild",
                    "sequence": sequence,
                    "smiles": smiles
                })

    # Generate CataPro input file
    catapro_input_df = pd.DataFrame(catapro_input)
    # Generate reverse mapping from SMILES to KEGG IDs as TSV 
    substrates_to_smiles_df = pd.DataFrame(list(substrates_to_smiles.items()), columns=['kegg_id', 'smiles'])

    return catapro_input_df, substrates_to_smiles_df


# --- Integrate CataPro predictions into kcat file ---


def integrate_catapro_predictions(kcat_df, substrates_to_smiles, catapro_predictions_df):
    """
    TODO: Write the documentation 
    """
    # Format the output to match the kcat_df
    # Convert pred_log10[kcat(s^-1)] to kcat(s^-1)
    catapro_predictions_df['kcat_s'] = 10 ** catapro_predictions_df['pred_log10[kcat(s^-1)]']
    catapro_predictions_df['uniprot_model'] = catapro_predictions_df['fasta_id'].str.replace('_wild', '', regex=False) # Extract UniProt ID
    # Match the SMILES to KEGG IDs using substrates_to_smiles
    catapro_predictions_df['substrates_kegg'] = catapro_predictions_df['smiles'].map(substrates_to_smiles.set_index('smiles')['kegg_id'])
    
    catapro_map = catapro_predictions_df.set_index(['uniprot_model', 'substrates_kegg'])['kcat_s'].to_dict()

    def get_min_pred_kcat(row):
        uniprot = row['uniprot_model']
        kegg_ids = str(row['substrates_kegg']).split(';')
        kcat_values = [
            catapro_map.get((uniprot, kegg_id))
            for kegg_id in kegg_ids
            if (uniprot, kegg_id) in catapro_map
        ]
        return min(kcat_values) if kcat_values else None  # If multiple substrates, take the minimum kcat value

    kcat_df['catapro_predicted_kcat_s'] = kcat_df.apply(get_min_pred_kcat, axis=1)
    return kcat_df


# --- Main ---


def run_catapro_part1(kcat_file_path, limit_matching_score, output_path, report=True):
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
    catapro_input_df, substrates_to_smiles_df = create_catapro_input_file(kcat_df, output_path)

    # Save the CataPro input file and substrates to SMILES mapping
    catapro_input_df.to_csv(output_path, sep=',', index=True)
    substrates_to_smiles_df.to_csv(output_path.replace('.csv', '_substrates_to_smiles.tsv'), sep='\t', index=False)


def run_catapro_part2(kcat_file_path, catapro_predictions_path, substrates_to_smiles_path, output_path, report=True):
    """
    TODO: Write the documentation
    """ 
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    substrates_to_smiles = pd.read_csv(substrates_to_smiles_path, sep='\t')
    catapro_predictions_df = pd.read_csv(catapro_predictions_path, sep=',')
    kcat_df = integrate_catapro_predictions(kcat_df, 
                                            substrates_to_smiles,
                                            catapro_predictions_df
                                            )
    
    # Save the output as a TSV file
    kcat_df.to_csv(output_path, sep='\t', index=False)
    # TODO: Generate report if required



if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))

    # Test : Integrate CataPro predictions into kcat file
    # kcat_df = pd.read_csv("output/ecoli_kcat_sabio.tsv", sep='\t')
    # substrates_to_smiles = pd.read_csv('in_progress/ml_test/substrates_to_smiles.tsv', sep='\t')
    # integrate_catapro_predictions(kcat_df, substrates_to_smiles, "in_progress/ml_test/catapro_output.csv", "in_progress/ml_test/ecoli_kcat_catapro.tsv")

    # Test : Main function
    run_catapro_part1("output/ecoli_kcat.tsv", 9, "in_progress/ml_test/catapro_input.csv")
    run_catapro_part2("output/ecoli_kcat_complete.tsv", 
                      "in_progress/ml_test/catapro_output.csv", 
                      "in_progress/ml_test/catapro_input_substrates_to_smiles.tsv", 
                      "output/ecoli_kcat_catapro.tsv")