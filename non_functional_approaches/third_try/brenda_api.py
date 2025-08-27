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

from non_functional_approaches.third_try.matching import exact_match_substrate, fuzzy_match_substrate, create_kcat_value
from wildkcat.utils.generate_reports import report_brenda


# TODO: Do a matching process not perfect but fonctional to begin with
# 
# 1 : match ec, substrate, enzyme, organism, variant, pH, temperature (SABIO-RK only)
# 2 : match ec, substrate, organism, variant, pH, temperature
# 3 : match ec, substrate, organism, variant
# 4 : match ec, substrate, organism
# 5 : match ec, substrate
# 6 : match ec


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
    if text is None or pd.isna(text):
        return None
    text = text.lower()
    if "wild" in text:  # wild-type, wildtype or wild type
        return "wildtype"
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


def find_best_match(brenda_df, kcat_dict, general_criterias):
    """
    TODO: Normalize the matching process across the different APIs.
    """
    ec_code = kcat_dict['ec_code']
    substrates_model = kcat_dict['substrates_name']

    # Filter by EC number
    candidates_base = brenda_df[brenda_df['ecNumber'] == ec_code].copy()
    if candidates_base.empty:
        return None, 8

    # Check criteria dynamically
    def check_general_criterias(row, keys_to_check):
        for key in keys_to_check:
            if key in general_criterias and not general_criterias[key](row[key]):
                return False
        return True

    # Substrate matching (exact, fallback fuzzy)
    def substrate_match(row_substrate):
        if exact_match_substrate(row_substrate, substrates_model):
            return True
        else:
            return False
        # fuzzy, score = fuzzy_match_substrate(row_substrate, substrates_model)
        # return bool(fuzzy)

    # Matching levels (strictest → loosest)
    matching_levels = [
        ['substrate', 'organism', 'variant', 'ph', 'temperature'],  # Level 2
        ['organism', 'variant', 'ph', 'temperature'],               # Level 3
        ['organism', 'variant'],                                    # Level 4
        ['organism'],                                               # Level 5
        []                                                          # Level 6
    ]

    for idx, criteria in enumerate(matching_levels, start=2):
        # Reset candidates each iteration to base EC filter
        candidates = candidates_base.copy()

        # Substrate filter
        if 'substrate' in criteria:
            candidates = candidates[candidates['substrate'].apply(substrate_match)]
            # Remove 'substrate' key since it's handled
            criteria = [c for c in criteria if c != 'substrate']

        # Other criteria
        if criteria:
            candidates = candidates[candidates.apply(lambda row: check_general_criterias(row, criteria), axis=1)]

        # Return best match if found
        if not candidates.empty:
            best_row = candidates.loc[candidates['turnoverNumber'].idxmax()]
            return best_row['turnoverNumber'], idx

    # No match found
    return None, 7


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

    if report:
        # Generate report
        report_brenda(kcat_df)


if __name__ == "__main__":
    # Test : Send a request to BRENDA API
    # df = get_turnover_number_brenda(
    #     ec_number="2.3.1.54",
    # )
    # df.to_csv("in_progress/brenda_test.tsv", sep='\t', index=False)

    # Test: Main function
    run_brenda('output/ecoli_kcat.tsv',
               'Escherichia coli',
               (20, 37),
               (6, 8)
        )