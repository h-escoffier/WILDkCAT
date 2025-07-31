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
        logging.warning('%s: No data found for the query.' % f"{ec_number}")
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


@lru_cache(maxsize=None)
def cached_get_turnover_number_brenda(ec_number):
    """Cached wrapper to avoid repeated BRENDA API calls."""
    return get_turnover_number_brenda(ec_number)


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


def extract_kcat_from_brenda(kcat_dict, general_criterias): 
    """
    TODO:
    """
    brenda_df = cached_get_turnover_number_brenda(kcat_dict['ec_code'])
    if brenda_df.empty: 
        return (None, 10) # No corresponding data for the EC code in BRENDA
    # Find the best match 
    return find_best_match(brenda_df, kcat_dict, general_criterias)
    

def run_brenda(kcat_file_path, organism, temperature_range, pH_range, variant = "wildtype", report=True):  # TODO: The run function could be shared with the other APIs
    """
    TODO: 
    """
    general_criterias = {
        "Organism": lambda x: x == organism,
        "Temperature": lambda x: temperature_range[0] <= x <= temperature_range[1],
        "pH": lambda x: pH_range[0] <= x <= pH_range[1],
        "Enzyme Variant": lambda x: x == variant,
    }

    # Read the kcat file
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')

    # Initialize new columns
    kcat_df['kcat'] = None
    kcat_df['matching_score'] = None

    for row in tqdm(kcat_df.itertuples(), total=len(kcat_df), desc="Processing BRENDA"):
        kcat_dict = row._asdict()

        # Extract kcat and matching score
        kcat, matching_score = extract_kcat_from_brenda(kcat_dict, general_criterias)

        # Assign results to the main dataframe
        kcat_df.loc[row.Index, 'kcat'] = kcat
        kcat_df.loc[row.Index, 'matching_score'] = matching_score

    kcat_df.to_csv("output/ecoli_kcat_brenda.tsv", sep='\t', index=False) # TODO: avoid to hardcode the output path
    logging.info("Output saved to 'output/ecoli_kcat_brenda.tsv'")


if __name__ == "__main__":
    # Test : Send a request to BRENDA API
    df = get_turnover_number_brenda(
        ec_number="2.7.1.11",
    )
    df.to_csv("in_progress/brenda_test.tsv", sep='\t', index=False)
