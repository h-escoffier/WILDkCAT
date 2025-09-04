import pandas as pd 
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from zeep import Client, Settings
from zeep.transports import Transport
from zeep.cache import InMemoryCache
from zeep.helpers import serialize_object
from dotenv import load_dotenv
from functools import lru_cache
import hashlib
import os 


from ..utils.manage_warnings import logger_retrieval as logger


load_dotenv()


# --- Setup --- 


def create_brenda_client(wsdl_url: str = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl") -> Client:
    """
    Creates and configures a persistent SOAP client for the BRENDA API.

    Parameters:
        wsdl_url (str): URL to the BRENDA WSDL file.

    Returns:
        zeep.Client: Configured SOAP client.
    """
    # Configure retry logic for network resilience
    session = Session()
    retry = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    # Set a custom User-Agent (some servers block default Python UA)
    session.headers.update({"User-Agent": "BRENDA-Client"})

    # Create zeep transport and settings
    transport = Transport(session=session, cache=InMemoryCache())
    settings = Settings(strict=False, xml_huge_tree=True) 

    return Client(wsdl_url, settings=settings, transport=transport)


def get_brenda_credentials() -> tuple[str, str]:
    """
    Retrieves and hashes BRENDA API credentials from environment variables.

    Returns:
        tuple[str, str]: (email, hashed_password)
    """
    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")

    if not email or not password:
        raise ValueError("BRENDA_EMAIL and BRENDA_PASSWORD environment variables must be set.")

    hashed_password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    return email, hashed_password


# --- BRENDA API ---


@lru_cache(maxsize=None)
def get_turnover_number_brenda(ec_number):
    """
    Queries the BRENDA SOAP API to retrieve turnover number values for a Enzyme Commission (EC) Number.

    Parameters:
        ec_number (str): EC number (e.g., '1.1.1.1').

    Returns:
        pd.DataFrame: A DataFrame containing turnover number entries.
    """

    email, hashed_password = get_brenda_credentials()
    client = create_brenda_client()

    # Define the parameters for the SOAP request

    parameters_kcat = [
        email,
        hashed_password,
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
        hashed_password,
        f'ecNumber*{ec_number}',
        "organism*",
        "sequenceCode*", 
        "commentary*", 
        "literature*",
        "textmining*"
    ]

    # print(client.service.__getattr__('getTurnoverNumber').__doc__)
    # print(client.service.__getattr__('getOrganism').__doc__)
    
    result_kcat = client.service.getTurnoverNumber(*parameters_kcat)
    result_organism = client.service.getOrganism(*parameters_org)
    
    # Format the response into a DataFrame
    data = serialize_object(result_kcat)
    data_organism = serialize_object(result_organism)

    if not data:
        logger.warning('%s: No data found for the query in BRENDA.' % f"{ec_number}")
        return pd.DataFrame()

    # Remove None values (-999)
    data = [entry for entry in data if entry.get('turnoverNumber') is not None and entry.get('turnoverNumber') != '-999']
    if data == []:
        logger.warning('%s: No valid data found for the query in BRENDA.' % f"{ec_number}")
        return pd.DataFrame()

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
    df["EnzymeVariant"] = df["commentary"].apply(get_variant)
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
    """
    Extracts the enzyme variant information from the commentary text.
    
    Parameters:
        text (str): Commentary text from BRENDA API response.

    Returns:
        str: The extracted enzyme variant information: wildtype, mutant, or None if not found.
    """
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
    """
    Queries the BRENDA SOAP API to retrieve cofactor information for a given Enzyme Commission (EC) number.

    Parameters:
        ec_number (str): EC number (e.g., '1.1.1.1').

    Returns:
        pd.DataFrame: A DataFrame containing turnover number entries.
    """
    # Call the SOAP API
    email, hashed_password = get_brenda_credentials()
    client = create_brenda_client()

    parameters_cofactor = [
        email,
        hashed_password,
        f'ecNumber*{ec_number}',
        "cofactor*", 
        "commentary*", 
        "organism*", 
        "ligandStructureId*", 
        "literature*"
    ]

    result_cofactor = client.service.getCofactor(*parameters_cofactor)
    data = serialize_object(result_cofactor)
    df = pd.DataFrame(data)
    if df.empty:
        return []
    cofactor = df['cofactor'].unique().tolist()
    return cofactor


# if __name__ == "__main__":
    # Test : Send a request to BRENDA API
    # df = get_turnover_number_brenda(ec_number="2.5.1.3")
    # df.to_csv("in_progress/api_output_test/brenda_test.tsv", sep='\t', index=False)

    # Test : Identify cofactor
    # df = get_cofactor("1.1.1.42")
    # df.to_csv("in_progress/api_output_test/brenda_cofactor.tsv", sep='\t', index=False)
