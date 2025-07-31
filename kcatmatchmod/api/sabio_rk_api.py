import requests
import logging
import pandas as pd 
from tqdm import tqdm
from io import StringIO
from functools import lru_cache


from kcatmatchmod.utils.matching import create_kcat_value
from kcatmatchmod.reports.generate_reports import report_sabio_rk


# TODO: Implement matching score calculation with a expert of the field


# ---------------------------------------------
# Retrieve turnover number (kcat) from SABIO-RK 
# for a given EC number & KEGG reaction ID
# ---------------------------------------------


def get_turnover_number_sabio(
        ec_number,
        kegg_reaction_id='',
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
    if kegg_reaction_id:
        query_parts.append(f'KeggReactionID:"{kegg_reaction_id}"')
    query_string = ' AND '.join(query_parts)
    query = {'format':'txt', 'q':query_string}

    # Make GET request
    request = requests.get(base_url, params=query)
    request.raise_for_status()  # Raise if 404 error
    if request.text == "no data found":
        logging.warning('%s: No data found for the query.' % f"{kegg_reaction_id} - {ec_number}")
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


@lru_cache(maxsize=None)
def cached_get_turnover_number_sabio(ec_number, kegg_reaction_id):
    """Cached wrapper to avoid repeated SABIO-RK API calls."""
    return get_turnover_number_sabio(ec_number, kegg_reaction_id)


# ---------------------------------------------
# Find substrates and products in SABIO-RK
# ---------------------------------------------


def to_set(s):
    return set(s.split(';')) if pd.notna(s) else set()


def find_rxn_direction(kcat_dict, sabio_df):
    """
    Determines the reaction direction in SABIO-RK data by matching substrate and product names.

    Parameters:
        kcat_dict (dict): Dictionary containing substrate and product names under 'substrates_name' and 'products_name'.
        sabio_df (pd.DataFrame): DataFrame containing SABIO-RK entries with 'Substrate' and 'Product' columns.

    Returns:
        pd.DataFrame: Filtered SABIO-RK DataFrame with matching substrates or products, or empty DataFrame if no match.
    """
    # Convert lists to sets
    substrates_kcat = to_set(kcat_dict['substrates_name'])
    products_kcat = to_set(kcat_dict['products_name'])

    # 1. Look for a substrate match
    sabio_df_substrate = sabio_df.copy()
    sabio_df_substrate['Substrate_set'] = sabio_df_substrate['Substrate'].apply(to_set)
    substrate_matches = sabio_df_substrate[
        sabio_df_substrate['Substrate_set'].apply(lambda s: not substrates_kcat.isdisjoint(s))
    ]
    if not substrate_matches.empty:
        # drop 'Substrate_set' column
        substrate_matches = substrate_matches.drop(columns=['Substrate_set'])
        return substrate_matches, None
    
    # 2. If no substrate match, look for a product match
    sabio_df_product = sabio_df.copy()
    sabio_df_product['Product_set'] = sabio_df_product['Product'].apply(to_set) 
    product_matches = sabio_df_product[
        sabio_df_product['Product_set'].apply(lambda p: not products_kcat.isdisjoint(p))
    ]
    if not product_matches.empty:
        product_matches = product_matches.drop(columns=['Substrate_set', 'Product_set'], errors='ignore')
        return product_matches, None

    # 3. If neither substrate nor product matches, check if the opposite direction matches
    sabio_df_substrate = sabio_df.copy()
    sabio_df_substrate['Substrate_set'] = sabio_df_substrate['Substrate'].apply(to_set)
    substrate_matches = sabio_df_substrate[
        sabio_df_substrate['Substrate_set'].apply(lambda s: not products_kcat.isdisjoint(s))
    ]
    if not substrate_matches.empty:
        return pd.DataFrame(), 8
    
    sabio_df_product = sabio_df.copy()
    sabio_df_product['Product_set'] = sabio_df_product['Product'].apply(to_set)
    product_matches = sabio_df_product[
        sabio_df_product['Product_set'].apply(lambda p: not substrates_kcat.isdisjoint(p))
    ]
    if not product_matches.empty:
        return pd.DataFrame(), 8

    # 4. If neither substrate nor product matches, log a warning
    kegg_rxn_id = sabio_df['KeggReactionID'].iloc[0] if 'KeggReactionID' in sabio_df.columns else 'Unknown'
    ec_code = sabio_df['ECNumber'].iloc[0] if 'ECNumber' in sabio_df.columns else 'Unknown'
    logging.warning('%s: No substrate or product matches between the two data sets.' % f"{kegg_rxn_id} - {ec_code}")
    return pd.DataFrame(), 9


# ---------------------------------------------
# Matching functions 
# ---------------------------------------------


def _apply_matching(kcat_dict, sabio_df, general_criterias, exclude_criteria_keys=None):
    """
    Generalized matching function to filter sabio_df based on kcat_dict and general_criterias,
    with optional exclusion of certain criteria keys.
    """
    mask = pd.Series(True, index=sabio_df.index)
    exclude_criteria_keys = exclude_criteria_keys or set()

    # Match UniProtKB_AC if available
    if 'uniprot_model' in kcat_dict and 'UniProtKB_AC' in sabio_df.columns:
        mask &= sabio_df['UniProtKB_AC'] == kcat_dict['uniprot_model']

    # Apply general criteria filtering
    for sabio_col, criteria_func in general_criterias.items():
        if sabio_col in exclude_criteria_keys:
            continue
        if sabio_col in sabio_df.columns:
            mask &= sabio_df[sabio_col].apply(criteria_func)

    return sabio_df[mask].copy()


def match_exact(kcat_dict, sabio_df, general_criterias):
    """Exact match using all criteria and UniProtKB_AC."""
    return _apply_matching(kcat_dict, sabio_df, general_criterias)


def match_relax_temp_pH(kcat_dict, sabio_df, general_criterias):
    """Match ignoring temperature and pH differences."""
    return _apply_matching(
        kcat_dict, sabio_df, general_criterias,
        exclude_criteria_keys={'Temperature', 'pH'}
    )


def match_relax_enzyme(kcat_dict, sabio_df, general_criterias):
    """Match ignoring enzyme (UniProtKB_AC) differences."""
    # Remove UniProtKB_AC matching by not applying it
    mask = pd.Series(True, index=sabio_df.index)
    exclude_criteria_keys = set()
    for sabio_col, criteria_func in general_criterias.items():
        if sabio_col in exclude_criteria_keys:
            continue
        if sabio_col in sabio_df.columns:
            mask &= sabio_df[sabio_col].apply(criteria_func)
    return sabio_df[mask].copy()


def match_relax_variant(kcat_dict, sabio_df, general_criterias):
    """Match ignoring variant (wildtype/mutant)."""
    return _apply_matching(
        kcat_dict, sabio_df, general_criterias,
        exclude_criteria_keys={'Enzyme Variant'}
    )


def match_ec_organism(kcat_dict, sabio_df, general_criterias):
    """Match EC code and organism only."""
    return _apply_matching(
        kcat_dict, sabio_df, general_criterias,
        exclude_criteria_keys={'Temperature', 'pH', 'Enzyme Variant'}
    )


def match_ec_only(kcat_dict, sabio_df, general_criterias):
    """Match only EC code (different organism, enzyme, temperature, pH and variant)."""
    mask = pd.Series(True, index=sabio_df.index)
    # No UniProtKB_AC, no organism, temperature, pH, variant
    return sabio_df[mask].copy()


# ---------------------------------------------
# Find best match
# ---------------------------------------------


def find_best_match(kcat_dict, sabio_df, general_criterias):
    """
    Attempt to find the best match in a hierarchical way:
    1. Perfect match
    2. Relax temperature/pH
    3. Relax enzyme
    4. Relax variant
    5. EC + organism
    6. EC only
    7. No match → return None (to be predicted later)
    """
    # 1. Perfect match
    match = match_exact(kcat_dict, sabio_df, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 1
    # 2. Relax temperature and pH
    match = match_relax_temp_pH(kcat_dict, sabio_df, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 2
    # 3. Relax the enzyme
    match = match_relax_enzyme(kcat_dict, sabio_df, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 3
    # 4. Relax enzyme variant
    match = match_relax_variant(kcat_dict, sabio_df, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 4
    # 5. Match EC code and organism
    match = match_ec_organism(kcat_dict, sabio_df, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 5
    # 6. Match EC code only
    match = match_ec_only(kcat_dict, sabio_df, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 6
    # 7. No match found
    logging.warning('%s: Matching failed' % (kcat_dict['KEGG_rxn_id'], kcat_dict['ec_code'])) # Should never happen
    return None, 7

        
# ---------------------------------------------
# Main functions
# ---------------------------------------------


def extract_kcat_from_sabio(kcat_dict, general_criterias):  
    """
    Extracts the kcat value from SABIO-RK for a given kcat_dict and general criteria.

    Parameters:
        kcat_dict (dict): Dictionary containing kcat information (EC code, KEGG reaction ID, substrates, products, etc.).
        general_criterias (dict): Dictionary of general criteria functions for filtering SABIO-RK entries.

    Returns:
        tuple: (kcat value or None, matching score)
    """
    sabio_df = cached_get_turnover_number_sabio(kcat_dict['ec_code'], kcat_dict['KEGG_rxn_id'])
    if sabio_df.empty: 
        return (None, 10) # No corresponding data for the EC code and KEGG reaction ID in SABIO-RK
    sabio_df, matching_score = find_rxn_direction(kcat_dict, sabio_df)
    if sabio_df.empty: 
        return (None, matching_score) # No match found in SABIO-RK for the given substrates/products
    # Find the best match
    return find_best_match(kcat_dict, sabio_df, general_criterias)


def run_sabio_rk(kcat_file_path, organism, temperature_range, pH_range, variant = "wildtype", report=True):
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

    for row in tqdm(kcat_df.itertuples(), total=len(kcat_df), desc="Processing SABIO-RK"):
        kcat_dict = row._asdict()

        # if 'KEGG_rxn_id' is missing, skip this row
        if 'KEGG_rxn_id' not in kcat_dict or pd.isna(kcat_dict['KEGG_rxn_id']):  # TODO: Should be process also 
            continue 

        # Extract kcat and matching score
        kcat, matching_score = extract_kcat_from_sabio(kcat_dict, general_criterias)
        
        # Assign results to the main dataframe
        kcat_df.loc[row.Index, 'kcat'] = kcat
        kcat_df.loc[row.Index, 'matching_score'] = matching_score

    kcat_df.to_csv("output/ecoli_kcat_sabio.tsv", sep='\t', index=False)
    logging.info("Output saved to 'output/ecoli_kcat_sabio.tsv'")

    if report:
        report_sabio_rk(kcat_df)


if __name__ == "__main__":
    # Test : Send a request to SABIO-RK API
    # df = get_turnover_number_sabio(
    #     ec_number="2.7.1.11",
    #     # kegg_reaction_id="R00209",
    # )
    # df.to_csv("in_progress/sabio_rk_test.tsv", sep='\t', index=False)

    # Test : Main function 
    run_sabio_rk('output/ecoli_kcat.tsv', 
                 'Escherichia coli',
                 (20, 37),
                 (6, 8)
                 )
