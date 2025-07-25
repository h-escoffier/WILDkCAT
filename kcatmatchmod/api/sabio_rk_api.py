import requests
import logging
import pandas as pd 
from io import StringIO
from functools import lru_cache
import pandas as pd
from tqdm import tqdm
import datetime
import os


# ---------------------------------------------
# Retrieve turnover number (kcat) from SABIO-RK 
# for a given EC number & KEGG reaction ID
# ---------------------------------------------


def get_turnover_number_sabio(
        ec_number,
        kegg_reaction_id: str = "",
        ):
    """
    TODO 
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
        logging.warning('%s: No data found for the query.' % query_string)
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
                                         'KineticMechanismType', 'Parameter']}

    # Make POST request
    request = requests.post(parameters, params=query, data=data_field)
    request.raise_for_status()

    # Format the response into a DataFrame
    df = pd.read_csv(StringIO(request.text), sep='\t')
    df = df[df['parameter.name'].str.lower() == 'kcat'].reset_index(drop=True) # Keep only kcat parameters
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
    TODO : I coded this to find the direction of the rxn bc I only have names :/// 
    """
    # Convert lists to sets
    substrates_kcat = to_set(kcat_dict['substrates_name'])
    products_kcat = to_set(kcat_dict['products_name'])

    # 1. Look for a substrate match
    sabio_df['Substrate_set'] = sabio_df['Substrate'].apply(to_set)
    substrate_matches = sabio_df[
        sabio_df['Substrate_set'].apply(lambda s: not substrates_kcat.isdisjoint(s))
    ]
    if not substrate_matches.empty:
        # drop 'Substrate_set' column
        substrate_matches = substrate_matches.drop(columns=['Substrate_set'])
        return substrate_matches

    # 2. If no substrate match, look for a product match
    sabio_df['Product_set'] = sabio_df['Product'].apply(to_set) 
    product_matches = sabio_df[
        sabio_df['Product_set'].apply(lambda p: not products_kcat.isdisjoint(p))
    ]
    if not product_matches.empty:
        product_matches = product_matches.drop(columns=['Substrate_set', 'Product_set'])
        return product_matches  # return all rows with product match

    # 3. If neither substrate nor product matches, log a warning
    kegg_rxn_id = sabio_df['KeggReactionID'].iloc[0] if 'KeggReactionID' in sabio_df.columns else 'Unknown'
    ec_code = sabio_df['ECNumber'].iloc[0] if 'ECNumber' in sabio_df.columns else 'Unknown'
    logging.warning('%s: No substrate or product matches between the two data sets.' % f"{kegg_rxn_id} - {ec_code}")
    return pd.DataFrame()


# ---------------------------------------------
# Helper functions for matching
# ---------------------------------------------


def _apply_matching(kcat_dict, sabio_df, mapping, general_criterias, exclude_mapping_keys=None, exclude_criteria_keys=None):
    """
    Generalized matching function to filter sabio_df based on mapping and general_criterias,
    with optional exclusion of certain keys.
    """
    mask = pd.Series(True, index=sabio_df.index)

    exclude_mapping_keys = exclude_mapping_keys or set()
    exclude_criteria_keys = exclude_criteria_keys or set()

    # Apply mapping-based filtering
    for kcat_key, sabio_col in mapping.items():
        if sabio_col in exclude_mapping_keys:
            continue
        if kcat_key not in kcat_dict or sabio_col not in sabio_df.columns:
            continue
        mask &= sabio_df[sabio_col] == kcat_dict[kcat_key]

    # Apply general criteria filtering
    for sabio_col, criteria_func in general_criterias.items():
        if sabio_col in exclude_criteria_keys:
            continue
        if sabio_col in sabio_df.columns:
            mask &= sabio_df[sabio_col].apply(criteria_func)

    return sabio_df[mask].copy()


def match_exact(kcat_dict, sabio_df, mapping, general_criterias):
    """Exact match using all mapping and criteria."""
    return _apply_matching(kcat_dict, sabio_df, mapping, general_criterias)


def match_relax_temp_pH(kcat_dict, sabio_df, mapping, general_criterias):
    """Match ignoring temperature and pH differences."""
    return _apply_matching(
        kcat_dict, sabio_df, mapping, general_criterias,
        exclude_criteria_keys={'Temperature', 'pH'}
    )


def match_relax_enzyme(kcat_dict, sabio_df, mapping, general_criterias):
    """Match ignoring enzyme differences."""
    return _apply_matching(
        kcat_dict, sabio_df, mapping, general_criterias,
        exclude_mapping_keys={'UniProtKB_AC'}
    )


def match_relax_variant(kcat_dict, sabio_df, mapping, general_criterias):
    """Match ignoring variant (wildtype/mutant)."""
    return _apply_matching(
        kcat_dict, sabio_df, mapping, general_criterias,
        exclude_criteria_keys={'Enzyme Variant'}
    )


def match_ec_organism(kcat_dict, sabio_df, mapping, general_criterias):
    """Match EC code and organism only."""
    return _apply_matching(
        kcat_dict, sabio_df, mapping, general_criterias,
        exclude_criteria_keys={'Temperature', 'pH', 'Enzyme Variant'}
    )


def match_ec_only(kcat_dict, sabio_df, mapping, general_criterias):
    """Match only EC code (different organism, enzyme, temperature, pH and variant)."""
    return _apply_matching(
        kcat_dict, sabio_df, mapping, general_criterias,
        exclude_mapping_keys={'UniProtKB_AC'},
        exclude_criteria_keys={'Organism', 'Temperature', 'pH', 'Enzyme Variant'}
    )


def create_kcat_value(df, method='mean'):
    """
    Aggregate kcat values from the DataFrame using the specified method.

    Parameters:
        df (pd.DataFrame): DataFrame containing SABIO-RK entries.
        method (str): Aggregation method: 'mean', 'max', or 'min'.

    Returns:
        float or None: Aggregated kcat value, or None if no valid values.
    """
    # Ensure the column exists and is numeric
    if 'parameter.startValue' not in df.columns:
        return None
    values = pd.to_numeric(df['parameter.startValue'], errors='coerce').dropna()
    if values.empty:
        return None
    if method == 'mean':
        return values.mean()
    elif method == 'max':
        return values.max()
    elif method == 'min':
        return values.min()
    else:
        raise ValueError("Invalid 'method' parameter. Choose from 'mean', 'max', or 'min'.")


# ---------------------------------------------
# Main function to find best match
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
    # Map kcat_dict keys to sabio_df columns
    key_mapping = {
            "ec_code": "ECNumber",
            "uniprot_model": "UniProtKB_AC",
            "substrates_name": "Substrates",
            "products_name": "Products",
        }
    # 1. Perfect match
    match = match_exact(kcat_dict, sabio_df, key_mapping, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 1
    # 2. Relax temperature and pH
    match = match_relax_temp_pH(kcat_dict, sabio_df, key_mapping, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 2
    # 3. Relax the enzyme
    match = match_relax_enzyme(kcat_dict, sabio_df, key_mapping, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 3
    # 4. Relax enzyme variant
    match = match_relax_variant(kcat_dict, sabio_df, key_mapping, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 4
    # 5. Match EC code and organism
    match = match_ec_organism(kcat_dict, sabio_df, key_mapping, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 5
    # 6. Match EC code only
    match = match_ec_only(kcat_dict, sabio_df, key_mapping, general_criterias)
    if not match.empty:
        kcat = create_kcat_value(match)
        return kcat, 6
    # 7. No match found
    logging.warning('%s: No match found in SABIO-RK - 7' % kcat_dict['ec_code']) # TODO: Should never appears
    return None, 7

        
# ---------------------------------------------
# Main functions
# ---------------------------------------------


def extract_kcat_from_sabio(kcat_dict, general_criterias):  
    """
    TODO
    """
    sabio_df = cached_get_turnover_number_sabio(kcat_dict['ec_code'], kcat_dict['KEGG_rxn_id'])
    if sabio_df.empty: 
        return (None, 9) # No corresponding data for the EC code and KEGG reaction ID in SABIO-RK
    sabio_df = find_rxn_direction(kcat_dict, sabio_df)
    if sabio_df.empty: 
        return (None, 8) # No match found in SABIO-RK for the given substrates/products
    # Find the best match
    return find_best_match(kcat_dict, sabio_df, general_criterias)


def run(kcat_file_path, organism, temperature_range, pH_range, variant = "wildtype", report=True):
    """
    TODO
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
        if 'KEGG_rxn_id' not in kcat_dict or pd.isna(kcat_dict['KEGG_rxn_id']):
            continue 

        # print(f"Processing {kcat_dict['ec_code']} with KEGG reaction ID {kcat_dict['KEGG_rxn_id']}...")

        # Extract kcat and matching score
        kcat, matching_score = extract_kcat_from_sabio(kcat_dict, general_criterias)
        
        # Assign results to the main dataframe
        kcat_df.loc[row.Index, 'kcat'] = kcat
        kcat_df.loc[row.Index, 'matching_score'] = matching_score

    kcat_df.to_csv("in_progress/kcat_with_scores.tsv", sep='\t', index=False)
    logging.info("Output saved to 'in_progress/kcat_with_scores.tsv'")
    if report:

        # Gather statistics
        total = len(kcat_df)
        matched = kcat_df['kcat'].notna().sum()
        match_percent = matched / total * 100 if total > 0 else 0
        score_counts = kcat_df['matching_score'].value_counts().sort_index()
        score_percent = (score_counts / total * 100).round(2)

        # Prepare HTML report
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html = f"""
        <html>
        <head><title>SABIO-RK kcat Matching Report</title></head>
        <body>
            <h1>SABIO-RK kcat Matching Report</h1>
            <p><b>Execution Time:</b> {now}</p>
            <p><b>Total entries:</b> {total}</p>
            <p><b>Matched kcat:</b> {matched} ({match_percent:.2f}%)</p>
            <h2>Matching Score Distribution</h2>
            <table border="1">
                <tr><th>Score</th><th>Count</th><th>Percent</th></tr>
        """
        for score, count in score_counts.items():
            percent = score_percent[score]
            html += f"<tr><td>{score}</td><td>{count}</td><td>{percent:.2f}%</td></tr>"
        html += """
            </table>
        </body>
        </html>
        """

        # Save report
        os.makedirs("reports", exist_ok=True)
        report_path = "reports/kcat_report.html"
        with open(report_path, "w") as f:
            f.write(html)
        logging.info(f"HTML report saved to '{report_path}'")

if __name__ == "__main__":
    # 1. Send a request to SABIO-RK API
    # df = get_turnover_number_sabio(
    #     ec_number="4.2.1.2",
    #     kegg_reaction_id="R01082",
    # )

    # Main function 
    run("output/ecoli_kcat.tsv", 
        'Escherichia coli',
        (20, 37),
        (6, 8), 
    )