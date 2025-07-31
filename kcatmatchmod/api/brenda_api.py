import logging
import pandas as pd 
from tqdm import tqdm
from zeep import Client, Settings
from zeep.helpers import serialize_object
from dotenv import load_dotenv
from functools import lru_cache
import hashlib
import os 
import re 


load_dotenv()


# ---------------------------------------------
# Retrieve turnover number (kcat) from BRENDA
# for a given EC number & KEGG reaction ID
# ---------------------------------------------


def get_turnover_number_brenda(
    ec_number: str,
):
    """
    Queries the BRENDA SOAP API to retrieve turnover number values for a given enzyme.

    Parameters:
        ec_number (str): Enzyme Commission (EC) number (e.g., '1.1.1.1').

    Returns:
        pd.DataFrame: A DataFrame containing turnover number entries.
    """

    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")

    # Call the SOAP API
    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    settings = Settings(strict=False)

    # print(client.service.__getattr__('getTurnoverNumber').__doc__)
    parameters = [
        email,
        password,
        f'ecNumber*{ec_number}',
        "turnoverNumber*", 
        "turnoverNumberMaximum*", 
        "substrate*", 
        "commentary*", 
        "organism*", 
        "ligandStructureId*", 
        "literature*"
    ]

    client = Client(wsdl, settings=settings)
    result = client.service.getTurnoverNumber(*parameters)
    

    # Format the response into a DataFrame
    data = serialize_object(result)
    if not data:
        return pd.DataFrame()

    # Remove None values (-999)
    data = [entry for entry in data if entry.get('turnoverNumber') is not None and entry.get('turnoverNumber') != '-999']
    df = pd.DataFrame(data)
    # Extract pH from commentary
    df["pH"] = df["commentary"].str.extract(r"pH\s*([\d\.]+)")
    # Extract temperature from commentary
    df["temperature"] = df["commentary"].str.extract(r"([\d\.]+)\?C")
    df["variant"] = df["commentary"].apply(get_variant)
    # Drop unnecessary columns
    df.drop(columns=["commentary", "ligandStructureId"], inplace=True, errors='ignore')
    return df


def get_variant(text):
    text = text.lower()
    if "wild" in text:  # wild-type ou wild type
        return "wild-type"
    elif any(word in text for word in ["mutant", "mutated", "mutation"]):
        return "mutant"
    return None


# ---------------------------------------------
# Matching functions 
# ---------------------------------------------


def matching(): 
    pass 


def match_organism_substrate_ph_temp(): 
    pass 


# ---------------------------------------------
# Find best match
# ---------------------------------------------


def find_best_match():
    pass 


# ---------------------------------------------
# Main functions
# ---------------------------------------------


def extract_kcat_from_brenda(): 
    pass 


def run_brenda(): 
    pass 


if __name__ == "__main__":
    # Test : Send a request to BRENDA API
    df = get_turnover_number_brenda(
        ec_number="2.7.1.11",
    )
    df.to_csv("in_progress/brenda_test.tsv", sep='\t', index=False)
