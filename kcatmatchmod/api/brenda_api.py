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
from kcatmatchmod.utils.generate_reports import report_api


load_dotenv()


# --- BRENDA API ---


@lru_cache(maxsize=None)
def get_turnover_number_brenda(ec_number):
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

    parameters_kcat = [
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

    parameters_org = [
        email,
        password,
        f'ecNumber*{ec_number}',
        "organism*",
        "sequenceCode*", 
        "commentary*", 
        "literature*",
        "textmining*"
    ]

    client = Client(wsdl, settings=settings)

    # print(client.service.__getattr__('getTurnoverNumber').__doc__)
    # print(client.service.__getattr__('getOrganism').__doc__)
    
    result_kcat = client.service.getTurnoverNumber(*parameters_kcat)
    result_organism = client.service.getOrganism(*parameters_org)
    
    # Format the response into a DataFrame
    data = serialize_object(result_kcat)
    data_organism = serialize_object(result_organism)

    if not data:
        logging.warning('%s: No data found for the query.' % f"{ec_number}")
        return pd.DataFrame()

    # Remove None values (-999)
    data = [entry for entry in data if entry.get('turnoverNumber') is not None and entry.get('turnoverNumber') != '-999']
    
    df = pd.DataFrame(data)
    df_org = pd.DataFrame(data_organism)

    # Format the organism response
    df_org.drop(columns=['commentary', 'textmining'], inplace=True, errors='ignore')
    
    # Merge on the literature column
    df_org['literature'] = df_org['literature'].apply(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x)
    df['literature'] = df['literature'].apply(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x)
    df = pd.merge(df, df_org, on=['literature', 'organism'], how='left')
    df.drop_duplicates(inplace=True)

    # Rename columns for consistency with other APIs
    df.rename(columns={
        'turnoverNumber': 'value',
        'sequenceCode' : 'UniProtKB_AC',
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
    df.drop(columns=["literature", "turnoverNumberMaximum", "parameter.endValue", "commentary", "ligandStructureId"], inplace=True, errors='ignore')
    
    # Remove the cofactor from the output 
    cofactor = get_cofactor(ec_number)
    # Drop the lines where the cofactor is not define
    df = df[~df['Substrate'].isin(cofactor)]   
    # Drop duplicates
    df.drop_duplicates(inplace=True)
    # Add a column for the db 
    df['db'] = 'brenda' 
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
def get_cofactor(ec_number):
    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")

    # Call the SOAP API
    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    settings = Settings(strict=False)

    parameters_cofactor = [
        email,
        password,
        f'ecNumber*{ec_number}',
        "cofactor*", 
        "commentary*", 
        "organism*", 
        "ligandStructureId*", 
        "literature*"
    ]

    client = Client(wsdl, settings=settings)
    result_cofactor = client.service.getCofactor(*parameters_cofactor)
    data = serialize_object(result_cofactor)
    df = pd.DataFrame(data)
    if df.empty:
        return []
    cofactor = df['cofactor'].unique().tolist()
    return cofactor


if __name__ == "__main__":
    # Test : Send a request to BRENDA API
    df = get_turnover_number_brenda(ec_number="2.7.1.11")
    df.to_csv("in_progress/api_output_test/brenda_test.tsv", sep='\t', index=False)

    # Test : Identify cofactor
    # df = get_cofactor("2.7.1.11")
    # df.to_csv("in_progress/api_output_test/brenda_cofactor.tsv", sep='\t', index=False)
