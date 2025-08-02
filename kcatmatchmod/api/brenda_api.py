import logging
import pandas as pd 
from tqdm import tqdm
from zeep import Client, Settings
from zeep.helpers import serialize_object
from dotenv import load_dotenv
from functools import lru_cache
import hashlib
import os 

from kcatmatchmod.utils.matching import find_best_match
from kcatmatchmod.reports.generate_reports import report_api


load_dotenv()


# --- BRENDA API ---


@lru_cache(maxsize=None)
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
    # Rename columns for consistency with other APIs
    df.rename(columns={
        'turnoverNumber': 'parameter.startValue',
        'turnoverNumberMaximum': 'parameter.endValue',
        'substrate': 'Substrate',
        'organism': 'Organism',
        'ecNumber': 'ECNumber'}, inplace=True) 

    # Extract pH from commentary
    df["pH"] = df["commentary"].str.extract(r"pH\s*([\d\.]+)")
    # Extract temperature from commentary
    df["Temperature"] = df["commentary"].str.extract(r"([\d\.]+)\?C")
    # Convert Temperature and pH to numeric, coercing errors to NaN
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['pH'] = pd.to_numeric(df['pH'], errors='coerce')
    # Extract enzyme variant from commentary
    df["Enzyme Variant"] = df["commentary"].apply(get_variant)
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


# --- Main ---


def extract_kcat_from_brenda(kcat_dict, general_criteria): 
    """
    TODO:
    """
    brenda_df = get_turnover_number_brenda(kcat_dict['ec_code'])
    if brenda_df.empty: 
        return (None, 10) # No corresponding data for the EC code in BRENDA
    # Find the best match 
    return find_best_match(kcat_dict, brenda_df, general_criteria, 'brenda')
    

def run_brenda(kcat_file_path, output_path, organism, temperature_range, pH_range, variant = "wildtype", report=True):  # TODO: The run function could be shared with the other APIs
    """
    TODO: 
    """
    general_criteria = {
        "Organism": organism,
        "Temperature": temperature_range,
        "pH": pH_range,
        "Enzyme Variant": variant,
    }

    # Read the kcat file
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')

    # Initialize new columns
    kcat_df['kcat'] = None
    kcat_df['matching_score'] = None

    for row in tqdm(kcat_df.itertuples(), total=len(kcat_df), desc="Processing BRENDA"):
        kcat_dict = row._asdict()

        # Extract kcat and matching score
        kcat, matching_score = extract_kcat_from_brenda(kcat_dict, general_criteria)

        # Assign results to the main dataframe
        kcat_df.loc[row.Index, 'kcat'] = kcat
        kcat_df.loc[row.Index, 'matching_score'] = matching_score

    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    if report:
        # Generate report
        report_api(kcat_df, "brenda")



if __name__ == "__main__":
    # Test : Send a request to BRENDA API
    # df = get_turnover_number_brenda(
    #     ec_number="2.7.1.11",
    # )
    # df.to_csv("in_progress/api_output_test/brenda_test.tsv", sep='\t', index=False)

    # Test: Main function
    # run_brenda('output/ecoli_kcat.tsv',
    #            'output/ecoli_kcat_brenda.tsv',
    #            'Escherichia coli',
    #            (20, 37),
    #            (6, 8)
    #     )
    
    # Test: Generate report
    df = pd.read_csv('output/ecoli_kcat_brenda.tsv', sep='\t')
    report_api(df, "brenda")