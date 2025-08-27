import logging
import pandas as pd 
from rapidfuzz import fuzz, process
import numpy as np


from wildkcat.utils.temperature import arrhenius_equation
from wildkcat.utils.organism import closest_org


# --- Substrates SABIO-RK --- 


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
        return substrate_matches
    
    # 2. If no substrate match, look for a product match
    sabio_df_product = sabio_df.copy()
    sabio_df_product['Product_set'] = sabio_df_product['Product'].apply(to_set) 
    product_matches = sabio_df_product[
        sabio_df_product['Product_set'].apply(lambda p: not products_kcat.isdisjoint(p))
    ]
    if not product_matches.empty:
        product_matches = product_matches.drop(columns=['Substrate_set', 'Product_set'], errors='ignore')
        return product_matches

    # 3. If neither substrate nor product matches, check if the opposite direction matches
    sabio_df_substrate = sabio_df.copy()
    sabio_df_substrate['Substrate_set'] = sabio_df_substrate['Substrate'].apply(to_set)
    substrate_matches = sabio_df_substrate[
        sabio_df_substrate['Substrate_set'].apply(lambda s: not products_kcat.isdisjoint(s))
    ]
    if not substrate_matches.empty:
        return pd.DataFrame()
    
    sabio_df_product = sabio_df.copy()
    sabio_df_product['Product_set'] = sabio_df_product['Product'].apply(to_set)
    product_matches = sabio_df_product[
        sabio_df_product['Product_set'].apply(lambda p: not substrates_kcat.isdisjoint(p))
    ]
    if not product_matches.empty:
        return pd.DataFrame()

    # 4. If neither substrate nor product matches, log a warning
    kegg_rxn_id = sabio_df['KeggReactionID'].iloc[0] if 'KeggReactionID' in sabio_df.columns else 'Unknown'
    ec_code = sabio_df['ECNumber'].iloc[0] if 'ECNumber' in sabio_df.columns else 'Unknown'
    logging.warning('%s: No substrate or product matches between the two data sets.' % f"{kegg_rxn_id} - {ec_code}")
    return pd.DataFrame()


    # sabio_df, matching_score = find_rxn_direction(kcat_dict, sabio_df)
    # if sabio_df.empty: 
    #     return (None, matching_score) # No match found in SABIO-RK for the given substrates/products
    # # Find the best match


# --- Check --- 
  

def check_enzyme(api_output, kcat_dict):
    """
    TODO: Write the documentation for this function.
    """
    enzyme = kcat_dict['uniprot_model']
    return api_output[api_output['UniProtKB_AC'] == enzyme]


def check_substrate(api_output, kcat_dict):
    """
    TODO: Write the documentation for this function.
    """
    results = []

    model_substrates = [s.strip().lower() for s in kcat_dict.get('substrates_name', '').split(';') if s.strip()]

    for api in api_output['db'].unique():
        db_subset = api_output[api_output['db'] == api]

        if api == 'sabio_rk':
            # 1. Try to match using KEGG reaction ID if provided
            if kcat_dict.get('KEGG_rxn_id'):
                kegg_rxn_matches = db_subset[db_subset['KeggReactionID'] == kcat_dict['KEGG_rxn_id']]
                if not kegg_rxn_matches.empty:
                    result = find_rxn_direction(kcat_dict, kegg_rxn_matches)
                    if not result.empty:
                        results.append(result)
                        continue  # If match found via KEGG ID, skip substrate matching
            
            # 2. Match using substrate names if no KEGG 
            def has_substrate_match(row_substrates):
                api_substrates = [s.strip().lower() for s in row_substrates.split(';') if s.strip()]
                return any(substrate in api_substrates for substrate in model_substrates)

            substrate_matches = db_subset[db_subset['Substrate'].apply(has_substrate_match)]
            if not substrate_matches.empty:
                results.append(substrate_matches)

        elif api == 'brenda':
            # Only match substrates
            substrate_matches = db_subset[db_subset['Substrate'].str.lower().isin(model_substrates)]
            if not substrate_matches.empty:
                results.append(substrate_matches)

        else:
            raise ValueError(f"Unknown API in 'db' column: {api}. Supported: 'sabio_rk', 'brenda'.")

    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return pd.DataFrame() 


def check_organism(api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    organism = general_criteria["Organism"]
    return api_output[api_output["Organism"] == organism]
    

def check_temp(api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    temp_min, temp_max = general_criteria["Temperature"]

    return api_output[
        (api_output["Temperature"].between(temp_min, temp_max))
    ]


def check_ph(api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    ph_min, ph_max = general_criteria["pH"]

    return api_output[
        (api_output["pH"].between(ph_min, ph_max))
    ]


def check_variant(api_output):
    """
    TODO: Write the documentation for this function.
    """
    return api_output[api_output["Enzyme Variant"] == 'wildtype']


# --- Matching --- 


def match_exact(kcat_dict, api_output, general_criteria):
    """
    TODO
    """
    df = api_output.copy()
    for check in [check_enzyme, check_substrate, check_organism, check_ph]:
        df = check(df, kcat_dict if check in [check_enzyme, check_substrate] else general_criteria)
        if df.empty:
            return None, None, None

    df_temp = check_temp(df, general_criteria)
    if not df_temp.empty:
        kcat, db = create_kcat_value(df_temp)
        return kcat, 1, db

    # If no temp match, try to estimate using Arrhenius equation
    kcat = arrhenius_equation(df, general_criteria)
    if kcat is not None:
        return kcat, 1.2, None
    return None, None, None


def relax_enzyme(kcat_dict, api_output, general_criteria):
    """
    TODO
    """
    df = api_output.copy()
    df = check_substrate(df, kcat_dict)
    if df.empty:
        return None, None, None

    df = check_ph(df, general_criteria)
    if df.empty:
        return None, None, None

    df = closest_org(kcat_dict, df)
    if not df.empty:
        kcat, db = create_kcat_value(df)
        return kcat, 2, db
    

    return None, None, None


def relax_organism(kcat_dict, api_output, general_criteria):
    """
    TODO
    """
    df = api_output.copy()
    df = check_substrate(df, kcat_dict)
    if df.empty:
        return None, None, None
    df = check_ph(df, general_criteria)
    if df.empty:
        return None, None, None

    df_temp = check_temp(df, general_criteria)
    if not df_temp.empty:
        kcat, db = create_kcat_value(df_temp)
        return kcat, 3, db

    kcat = arrhenius_equation(df, general_criteria)
    if kcat is not None:
        return kcat, 3.2, None
    return None, None, None


def relax_temp_and_ph(kcat_dict, api_output, general_criteria):
    """
    Relaxes temperature and pH constraints: matches enzyme, substrate, and organism only.
    """
    df = api_output.copy()
    for check in [check_enzyme, check_substrate, check_organism]:
        df = check(df, kcat_dict if check == check_enzyme or check == check_substrate else general_criteria)
        if df.empty:
            return None, None, None

    kcat, db = create_kcat_value(df)
    return kcat, 4, db


def just_ec(kcat_dict, api_output, general_criteria):
    """
    Only matches EC number, returns any available kcat value.
    """
    kcat, db = create_kcat_value(api_output)
    return kcat, 5, db


# --- Organism ---


def find_closest_organism(): 
    pass


# --- Iteration --- 


def find_best_match(kcat_dict, api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    # Keep only wildtype entries
    api_output = check_variant(api_output)

    # 1. Perfect match (only for SABIO-RK)
    kcat, matching_value, db = match_exact(kcat_dict, api_output, general_criteria)
    if kcat is not None:
        return kcat, matching_value, db

    # 2. Relax the organism
    kcat, matching_value, db = relax_organism(kcat_dict, api_output, general_criteria)
    if kcat is not None:
        return kcat, matching_value, db
    
    # 3. Relax the enzyme
    kcat, matching_value, db = relax_enzyme(kcat_dict, api_output, general_criteria)
    if kcat is not None:
        return kcat, matching_value, db

    # 4. Relax temperature and pH
    kcat, matching_value, db = relax_temp_and_ph(kcat_dict, api_output, general_criteria)
    if kcat is not None:
        return kcat, matching_value, db

    #5. Just EC
    kcat, matching_value, db = just_ec(kcat_dict, api_output, general_criteria)
    if kcat is not None:
        return kcat, matching_value, db

    return None, 9, None


# --- Generate kcat value --- 


def create_kcat_value(df, method='min'):
    """
    Aggregate kcat values from the DataFrame using the specified method and return the value along with its source database.

    Parameters:
        df (pd.DataFrame): DataFrame containing SABIO-RK or BRENDA entries.
        method (str): Aggregation method: 'max', or 'min'.

    Returns:
        tuple: (aggregated kcat value (float), source database (str)), or (None, None) if no valid values.
    """
    if 'value' not in df.columns or 'db' not in df.columns:
        return None, None
    values = pd.to_numeric(df['value'], errors='coerce')
    valid = df.loc[values.notna()]
    if valid.empty:
        return None, None
    if method == 'max':
        idx = values.idxmax()
    elif method == 'min':
        idx = values.idxmin()
    else:
        raise ValueError("Invalid 'method' parameter. Choose from 'max', or 'min'.")
    kcat_value = values.loc[idx]
    db_source = df.loc[idx, 'db']
    return kcat_value, db_source



# if __name__ == "__main__":
#     # Test : Test the matching function 
#     kcat_dict = {'Index': 0, 
#                  'rxn': 'PFK', 'KEGG_rxn_id': None, 
#                  'ec_code': '2.7.1.11', 
#                  'direction': 'forward', 
#                  'substrates_name': 'ATP C10H12N5O13P3;D-Fructose 6-phosphate', 
#                  'substrates_kegg': 'C00002;C05345', 'products_name': 'ADP C10H12N5O10P2;D-Fructose 1,6-bisphosphate;H+', 
#                  'products_kegg': 'C00008;C00354;C00080', 
#                  'genes_model': 'b3916', 'uniprot_model': 'P0A796', 
#                  'kegg_genes': 'b1723;b3916', 'intersection_genes': 'b3916', 
#                  'kcat': None, 'matching_score': None}
    
#     general_criteria = {
#         "Organism": 'Escherichia coli',
#         "Temperature": (20, 37),
#         "pH": (6, 8),
#         "Enzyme Variant": 'wildtype',
#     }