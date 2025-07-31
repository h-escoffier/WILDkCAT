import requests
import pandas as pd 
import re
from functools import lru_cache
import logging


# ---------------------------------------------
# Retrieve turnover number (kcat) from UniProt
# for a given UniProt ID
# ---------------------------------------------


def get_turnover_number_uniprot(uniprot_id):
    """
    Queries the UniProt REST API to retrieve turnover number values for a given UniProt ID.

    Parameters:
        uniprot_id (str): UniProt accession ID (e.g., 'P37689').

    Returns:
        pd.DataFrame: A DataFrame containing turnover number entries.
    """
    params = {
        "fields": [
            "kinetics"
        ]
    }

    headers = {
        "accept": "application/json"
    }

    base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    response = requests.get(base_url, headers=headers, params=params)
    
    if not response.ok:
        response.raise_for_status()
        logging.warning(f"Failed to retrieve data from UniProt API: {response.status_code} - {uniprot_id}")
        return pd.DataFrame()  # Return empty DataFrame on error

    data = response.json()

    kcat_entries = []

    # Extract kcat values from the response 
    comments = data.get("comments", [])
    for comment in comments:
        kinetic_params = comment.get("kineticParameters", {})
        note = kinetic_params.get("note", {})
        texts = note.get("texts", [])

        for text_block in texts:
            value = text_block.get("value", "")
            
            # Extract kcat and substrate
            matches = re.findall(r'kcat is (\d+(?:\.\d+)?) sec\(-1\) .*? with ([\w-]+) as substrate', value)
            for match in matches:
                kcat_value = float(match[0])
                substrate = match[1]

                # Extract pH
                ph_match = re.search(r'at pH (\d+\.?\d*)', value)
                ph_value = float(ph_match.group(1)) if ph_match else None

                # Extract temperature
                temp_match = re.search(r'(\d+) degrees Celsius', value)
                temp_value = float(temp_match.group(1)) if temp_match else None

                kcat_entries.append({
                    "UniProt_ID": uniprot_id,
                    "substrate": substrate,
                    "kcat": kcat_value,
                    "pH": ph_value,
                    "temperature": temp_value
                })

    uniprot_df = pd.DataFrame(kcat_entries)
    if uniprot_df.empty:
        logging.warning(f"No kcat data found for UniProt ID: {uniprot_id}")
        return pd.DataFrame()

    return uniprot_df


@lru_cache(maxsize=None)
def cached_get_turnover_number_uniprot(uniprot_id):
    return get_turnover_number_uniprot(uniprot_id)


# ---------------------------------------------
# Matching functions 
# ---------------------------------------------

def matching(): 
    pass 


# ---------------------------------------------
# Find best match
# ---------------------------------------------


def find_best_match():
    pass


# ---------------------------------------------
# Main functions 
# ---------------------------------------------


def run_uniprot(): 
    pass 


if __name__ == "__main__":
    # Test : Send a request to UniProt API
    uniprot_id = "P25437"
    df = get_turnover_number_uniprot(uniprot_id)
    df.to_csv("in_progress/uniprot_test.tsv", sep='\t', index=False)