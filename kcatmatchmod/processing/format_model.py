import re 
import logging
import requests
import pandas as pd 
from tqdm import tqdm
from functools import lru_cache

from kcatmatchmod.utils.model_function import read_model
# TODO: from kcatmatchmod.reports.generate_reports import report_format_model


# I don't know if this make sense ? 
# Find the KEGG rxn ID 
# Identify the EC-Code and check if the ec code is valid 
# If an EC code is not valid, find the new one 
# Check the GPR rules and check that the genes are matching the one in KEGG 
# Update the model with the new EC code and the genes 


def find_new_rxn_id(rxn_id): 
    """
    TODO: 
    """
    pass 


def find_new_ec_code(ec_code): 
    """
    TODO: 
    """
    url = f'https://rest.kegg.jp/list/{ec_code}'
    response = requests.get(url)

    if response.status_code != 200:
        logging.error(f"Failed to retrieve data from KEGG API: {response.status_code}")
        return None

    if "Transferred to" in response.text:
        match = re.search(r'Transferred to (.+)', response.text)
        if match:
            ec_codes = re.findall(r'\d+\.\d+\.\d+\.\d+', match.group(1))
            return ec_codes
    return []
    

def update_ec_code(): 
    pass


def find_inconsistent_gpr_rules(): 
    pass


def run_update_model(): 
    pass 


if __name__ == "__main__":
    pass

