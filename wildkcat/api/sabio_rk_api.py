import requests
import pandas as pd 
from io import StringIO
from functools import lru_cache 

from ..utils.manage_warnings import logger_retrieval as logger


# --- Sabio-RK API ---


@lru_cache(maxsize=None)
def get_turnover_number_sabio(ec_number):
    """
    Retrieve turnover number (kcat) data from SABIO-RK for a given EC number and KEGG reaction ID.

    Parameters:
        ec_number (str): Enzyme Commission number.
        kegg_reaction_id (str): KEGG reaction ID.

    Returns:
        pd.DataFrame: DataFrame containing SABIO-RK entries for kcat.
    """
    base_url = 'https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'
    parameters = 'https://sabiork.h-its.org/entry/exportToExcelCustomizable'
    entryIDs = []

    # -- Retrieve entryIDs --
    query = {'format': 'txt', 'q': f'Parametertype:"kcat" AND ECNumber:"{ec_number}"'}
    
    # Make GET request
    request = requests.get(base_url, params=query)
    request.raise_for_status()
    if request.text == "no data found":
        logger.warning('%s: No data found for the query in SABIO-RK.' % f"{ec_number}")
        return pd.DataFrame()  # Return empty DataFrame if no data found

    entryIDs = [int(x) for x in request.text.strip().split('\n')]
    # Retrieve informations matching the entryIDs
    data_field = {'entryIDs[]': entryIDs}
    # Possible fields to retrieve:
    # EntryID, Reaction, Buffer, ECNumber, CellularLocation, UniProtKB_AC, Tissue, Enzyme Variant, Enzymename, Organism
    # Temperature, pH, Activator, Cofactor, Inhibitor, KeggReactionID, KineticMechanismType, Other Modifier, Parameter,
    # Pathway, Product, PubMedID, Publication, Rate Equation, SabioReactionID, Substrate
    query = {'format':'tsv', 'fields[]':['EntryID', 'ECNumber', 'KeggReactionID', 'Reaction', 'Substrate', 'Product', 
                                         'UniProtKB_AC', 'Organism', 'Enzyme Variant', 'Temperature', 'pH', 
                                         'Parameter']}

    # Make POST request
    request = requests.post(parameters, params=query, data=data_field)
    request.raise_for_status()

    # Format the response into a DataFrame
    df = pd.read_csv(StringIO(request.text), sep='\t')
    df = df[df['parameter.name'].str.lower() == 'kcat'].reset_index(drop=True) # Keep only kcat parameters
    # Convert Temperature and pH to numeric, coercing errors to NaN
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['pH'] = pd.to_numeric(df['pH'], errors='coerce')
    # Drop unnecessary columns
    df.drop(columns=['EntryID', 'parameter.name', 'parameter.type', 'parameter.associatedSpecies', 
                     'parameter.endValue', 'parameter.standardDeviation'], inplace=True, errors='ignore')
    # Drop duplicates based on normalized Substrate and Product sets
    df["Substrate_set"] = df["Substrate"].fillna("").str.split(";").apply(lambda x: tuple(sorted(s.strip() for s in x if s.strip())))
    df["Product_set"] = df["Product"].fillna("").str.split(";").apply(lambda x: tuple(sorted(s.strip() for s in x if s.strip())))
    dedup_cols = [col for col in df.columns if col not in ["Substrate", "Product"]]
    df = df.drop_duplicates(subset=dedup_cols + ["Substrate_set", "Product_set"], keep="first")
    df = df.drop(columns=["Substrate_set", "Product_set"])
    # Rename columns for consistency
    df.rename(columns={
        'ECNumber': 'ECNumber',
        'KeggReactionID': 'KeggReactionID',
        'Reaction': 'Reaction',
        'Substrate': 'Substrate',
        'Product': 'Product',
        'UniProtKB_AC': 'UniProtKB_AC',
        'Organism': 'Organism',
        'Enzyme Variant': 'EnzymeVariant',
        'Temperature': 'Temperature',
        'pH': 'pH',
        'parameter.startValue': 'value',
        'parameter.unit': 'unit'
    }, inplace=True)
    # Add a column for the db
    df['db'] = 'sabio_rk'
    return df

        
# if __name__ == "__main__":
    # Test : Send a request to SABIO-RK API
    # df = get_turnover_number_sabio(ec_number="1.1.1.42")
    # df.to_csv("in_progress/api_output_test/sabio_rk_test.tsv", sep='\t', index=False)