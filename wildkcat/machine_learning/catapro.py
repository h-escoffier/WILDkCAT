import logging
import pandas as pd 
from tqdm import tqdm
import re
from functools import lru_cache

from wildkcat.api.api_utilities import safe_requests_get, retry_api
from wildkcat.api.uniprot_api import convert_uniprot_to_sequence, identify_catalytic_enzyme
from wildkcat.api.brenda_api import get_cofactor


# TODO: Add a list of cofactors 


# --- API ---


def convert_kegg_compound_to_sid(kegg_compound_id) -> str | None:
    """
    Convert the KEGG compound ID to the PubChem Substance ID (SID).

    Parameters:
        kegg_compound_id (str): KEGG compound ID.

    Returns:
        str: The PubChem SID if found, otherwise None.
    """
    url = f"https://rest.kegg.jp/conv/pubchem/compound:{kegg_compound_id}"
    safe_get_with_retry = retry_api(max_retries=2, backoff_factor=2)(safe_requests_get)
    response = safe_get_with_retry(url)

    if response is None:
        return None

    if response.status_code != 200:
        return None

    match = re.search(r'pubchem:\s*(\d+)', response.text)
    sid = match.group(1) if match else None
    return sid


def convert_sid_to_cid(sid) -> int | None:
    """
    Converts a PubChem Substance ID (SID) to the corresponding Compound ID (CID).

    Parameters:
        sid (str): PubChem Substance ID.

    Returns:
        int or None: The corresponding PubChem Compound ID (CID), or None if not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/JSON"
    safe_get_with_retry = retry_api(max_retries=1, backoff_factor=2)(safe_requests_get)
    response = safe_get_with_retry(url)
    
    if response is None:
        return None

    if response.status_code == 200:
        try:
            cid = response.json()['InformationList']['Information'][0]['CID'][0]
        except (KeyError, IndexError):
            cid = None
    return cid


def convert_cid_to_smiles(cid) -> list | None:    
    """
    Converts a PubChem Compound ID (CID) to its corresponding SMILES representation.

    Parameters:
        cid (str): PubChem Compound ID.

    Returns:
       list or None: A list of SMILES strings if found, otherwise None.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/smiles/txt"
    try:
        safe_get_with_retry = retry_api(max_retries=1, backoff_factor=2)(safe_requests_get)
        response = safe_get_with_retry(url)

        if response is None:
            return None

        response.raise_for_status()
        smiles = response.text.strip().split('\n')
        return smiles
    except:
        return None


@lru_cache(maxsize=None)
def convert_kegg_to_smiles(kegg_compound_id) -> list | None:
    """
    Convert the KEGG compound ID to the PubChem Compound ID (CID).

    Parameters:
        kegg_compound_id (str): KEGG compound ID.

    Returns:
        list or None: A list of SMILES strings if found, otherwise None.
    """
    sid = convert_kegg_compound_to_sid(kegg_compound_id)
    if sid is None:
        logging.warning('%s: Failed to retrieve SID for KEGG compound ID' % (kegg_compound_id))
        return None
    cid = convert_sid_to_cid(sid)
    if cid is None:
        logging.warning('%s: Failed to retrieve CID for KEGG compound ID' % (kegg_compound_id))
        return None
    smiles = convert_cid_to_smiles(cid)
    if smiles is None:
        logging.warning('%s: Failed to retrieve SMILES for KEGG compound ID' % (kegg_compound_id))
        return None
    return smiles
    

# --- Create CataPro input file ---


def create_catapro_input_file(kcat_df):
    """
    Generate CataPro input file and a mapping of substrate KEGG IDs to SMILES.

    Parameters: 
        kcat_df (pd.DataFrame): Input DataFrame containing kcat information.

    Returns:
        catapro_input_df (pd.DataFrame): DataFrame for CataPro input.
        substrates_to_smiles (dict): Mapping KEGG ID <-> SMILES.
    """
    catapro_input = []
    substrates_to_smiles = {}

    counter_no_catalytic, counter_kegg_no_matching, counter_rxn_covered, counter_cofactor = 0, 0, 0, 0
    for _, row in tqdm(kcat_df.iterrows(), total=len(kcat_df), desc="Generating CataPro input"):
        uniprot = row['uniprot']
        ec_code = row['ec_code']

        if len(uniprot.split(';')) > 1:       
            catalytic_enzyme = identify_catalytic_enzyme(uniprot, ec_code)
            if catalytic_enzyme is None or (";" in str(catalytic_enzyme)):
                counter_no_catalytic += 1
                continue
            else: 
                uniprot = catalytic_enzyme
                
        # If the number of KEGG Compound IDs is not matching the number of names, continue 
        if len([s for s in row['substrates_kegg'].split(';') if s]) != len(row['substrates_name'].split(';')):
            logging.warning(f"Number of KEGG compounds IDs does not match number of names for {ec_code}: {uniprot}.")
            counter_kegg_no_matching += 1
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
            if name.lower() in [c.lower() for c in cofactor]:  # TODO: Should we add a warning if no cofactor is found for a reaction? 
                counter_cofactor += 1
                continue
            smiles = convert_kegg_to_smiles(kegg_compound_id)
            if smiles is not None:
                smiles_str = smiles[0]  # TODO: If multiple SMILES, take the first one ? 
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
        
        counter_rxn_covered += 1

    # Generate CataPro input file
    catapro_input_df = pd.DataFrame(catapro_input)
    # Remove duplicates
    before_duplicates_filter = len(catapro_input_df)
    catapro_input_df = catapro_input_df.drop_duplicates().reset_index(drop=True)
    nb_lines_dropped = before_duplicates_filter - len(catapro_input_df)
    # Generate reverse mapping from SMILES to KEGG IDs as TSV
    substrates_to_smiles_df = pd.DataFrame(list(substrates_to_smiles.items()), columns=['kegg_id', 'smiles'])

    report_statistics = {
        "rxn_covered": counter_rxn_covered,
        "cofactor_identified": counter_cofactor,
        "no_catalytic": counter_no_catalytic,
        "kegg_no_matching": counter_kegg_no_matching,
        "duplicates_enzyme_substrates": nb_lines_dropped,
    }

    return catapro_input_df, substrates_to_smiles_df, report_statistics


# --- Integrate CataPro predictions into kcat file ---


def integrate_catapro_predictions(kcat_df, substrates_to_smiles, catapro_predictions_df) -> pd.DataFrame:
    """
    Integrates Catapro predictions into an kcat file.
    If multiple values are provided for a single combination of EC, Enzyme, Substrate, the minimum value is taken.

    Parameters:
        kcat_df (pd.DataFrame): Input DataFrame containing kcat information.
        substrates_to_smiles (pd.DataFrame): DataFrame mapping KEGG ID <-> SMILES.
        catapro_predictions_df (pd.DataFrame): DataFrame containing Catapro model predictions

    Returns:
        pd.DataFrame: The input kcat_df with an additional column 'catapro_predicted_kcat_s' containing
            the integrated Catapro predicted kcat(s^-1) values.
    """
    # Convert pred_log10[kcat(s^-1)] to kcat(s^-1)
    catapro_predictions_df['kcat_s'] = 10 ** catapro_predictions_df['pred_log10[kcat(s^-1)]']
    catapro_predictions_df['uniprot'] = catapro_predictions_df['fasta_id'].str.replace('_wild', '', regex=False) # Extract UniProt ID
    
    # Match the SMILES to KEGG IDs using substrates_to_smiles
    # If multiple KEGG IDs are found for a single SMILES, they are concatenated
    smiles_to_kegg = (
        substrates_to_smiles.groupby('smiles')['kegg_id']
        .apply(lambda x: ';'.join(sorted(set(x))))
    )
    catapro_predictions_df['substrates_kegg'] = catapro_predictions_df['smiles'].map(smiles_to_kegg)
    
    catapro_map = catapro_predictions_df.set_index(['uniprot', 'substrates_kegg'])['kcat_s'].to_dict()

    def get_min_pred_kcat(row):
        uniprot = row['uniprot']
        kegg_ids = str(row['substrates_kegg']).split(';')
        kcat_values = [
            catapro_map.get((uniprot, kegg_id))
            for kegg_id in kegg_ids
            if (uniprot, kegg_id) in catapro_map
        ]
        return min(kcat_values) if kcat_values else None  # If multiple substrates, take the minimum kcat value

    kcat_df['catapro_predicted_kcat_s'] = kcat_df.apply(get_min_pred_kcat, axis=1)
    return kcat_df


# if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))

    # Test : Integrate CataPro predictions into kcat file
    # kcat_df = pd.read_csv("output/ecoli_kcat_sabio.tsv", sep='\t')
    # substrates_to_smiles = pd.read_csv('in_progress/ml_test/substrates_to_smiles.tsv', sep='\t')
    # integrate_catapro_predictions(kcat_df, substrates_to_smiles, "in_progress/ml_test/catapro_output.csv", "in_progress/ml_test/ecoli_kcat_catapro.tsv")