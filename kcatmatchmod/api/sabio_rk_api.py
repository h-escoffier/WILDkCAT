import requests
import logging
import pandas as pd 
from tqdm import tqdm
from io import StringIO
from functools import lru_cache

from kcatmatchmod.api.api_utilities import retry_api
from kcatmatchmod.utils.matching import find_best_match
from kcatmatchmod.reports.generate_reports import report_api


# TODO: Implement matching score calculation with a expert of the field


# --- Sabio-RK API ---


@lru_cache(maxsize=None)
@retry_api(max_retries=4, backoff_factor=2)
def get_turnover_number_sabio(
        ec_number
        ):
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
    query_parts = ['Parametertype:"kcat"']
    if ec_number:
        query_parts.append(f'ECNumber:"{ec_number}"')
    query_string = ' AND '.join(query_parts)
    query = {'format':'txt', 'q':query_string}

    # Make GET request
    request = requests.get(base_url, params=query)
    request.raise_for_status()  # Raise if 404 error
    if request.text == "no data found":
        logging.warning('%s: No data found for the query.' % f"{ec_number}")
        return pd.DataFrame()  # Return empty DataFrame if no data found

    entryIDs = [int(x) for x in request.text.strip().split('\n')]

    # -- Retrieve informations matching the entryIDs -- 
    data_field = {'entryIDs[]': entryIDs}
    # Possible fields to retrieve:
    # EntryID, Reaction, Buffer, ECNumber, CellularLocation, UniProtKB_AC, Tissue, Enzyme Variant, Enzymename, Organism
    # Temperature, pH, Activator, Cofactor, Inhibitor, KeggReactionID, KineticMechanismType, Other Modifier, Parameter,
    # Pathway, Product, PubMedID, Publication, Rate Equation, SabioReactionID, Substrate
    query = {'format':'tsv', 'fields[]':['EntryID', 'ECNumber', 'KeggReactionID', 'Reaction', 'Substrate', 'Product', 
                                         'UniProtKB_AC', 'Organism', 'Enzyme Variant', 'Temperature', 'pH', 
                                         'Activator', 'Cofactor', 'Inhibitor',
                                         'KineticMechanismType', 'Parameter']}

    # Make POST request
    request = requests.post(parameters, params=query, data=data_field)
    request.raise_for_status()

    # Format the response into a DataFrame
    df = pd.read_csv(StringIO(request.text), sep='\t')
    df = df[df['parameter.name'].str.lower() == 'kcat'].reset_index(drop=True) # Keep only kcat parameters
    # Convert Temperature and pH to numeric, coercing errors to NaN
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['pH'] = pd.to_numeric(df['pH'], errors='coerce')
    return df

        
# --- Main ---


def extract_kcat_from_sabio(kcat_dict, general_criteria):  
    """
    Extracts the kcat value from SABIO-RK for a given kcat_dict and general criteria.

    Parameters:
        kcat_dict (dict): Dictionary containing kcat information (EC code, KEGG reaction ID, substrates, products, etc.).
        general_criteria (dict): Dictionary of general criteria functions for filtering SABIO-RK entries.

    Returns:
        tuple: (kcat value or None, matching score)
    """
    sabio_df = get_turnover_number_sabio(kcat_dict['ec_code'])
    if sabio_df.empty: 
        return None, 10 # No corresponding data for the EC code in SABIO-RK
    return find_best_match(kcat_dict, sabio_df, general_criteria, 'sabio_rk')


def run_sabio_rk(kcat_file_path, output_path, organism, temperature_range, pH_range, variant = "wildtype", report=True):
    """
    Runs the SABIO-RK kcat matching pipeline for a given input file and criteria.

    Parameters:
        kcat_file_path (str): Path to the input TSV file containing kcat data.
        organism (str): Organism name to match in SABIO-RK entries.
        temperature_range (tuple): Acceptable temperature range (min, max) for matching.
        pH_range (tuple): Acceptable pH range (min, max) for matching.
        variant (str, optional): Enzyme variant to match (default is "wildtype").
        report (bool, optional): If True, generates an HTML report of matching statistics.

    Returns:
        None. Outputs a TSV file with matched kcat values and an optional HTML report.
    """
    # Create a dict with the general criterias
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

    # Retrieve kcat values from SABIO-RK
    for row in tqdm(kcat_df.itertuples(), total=len(kcat_df), desc="Processing SABIO-RK"):
        kcat_dict = row._asdict()

        # Extract kcat and matching score
        kcat, matching_score = extract_kcat_from_sabio(kcat_dict, general_criteria)
        
        # Assign results to the main dataframe
        kcat_df.loc[row.Index, 'kcat'] = kcat
        kcat_df.loc[row.Index, 'matching_score'] = matching_score

    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    if report:
        report_api(kcat_df, 'sabio_rk')
    
    return kcat_df


if __name__ == "__main__":
    # Test : Send a request to SABIO-RK API
    # df = get_turnover_number_sabio(
    #     ec_number="2.7.1.11",
    #     # kegg_reaction_id="R00209",
    # )
    # df.to_csv("in_progress/api_output_test/sabio_rk_test.tsv", sep='\t', index=False)

    # Test : Main function 
    # run_sabio_rk('output/ecoli_kcat.tsv', 
    #              'output/ecoli_kcat_sabio.tsv',
    #              'Escherichia coli',
    #              (20, 37),
    #              (6, 8)
    #              )

    # Test: Generate report
    df = pd.read_csv('output/ecoli_kcat_sabio.tsv', sep='\t')
    report_api(df, 'sabio_rk')