import logging
import pandas as pd 
from rapidfuzz import fuzz, process



# ---------------------------------------------
# Matching functions - SABIO-RK specific
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


# ---------------------------------------------
# Matching substrates
# ---------------------------------------------



def exact_match_substrate(substrate_api, substrates_model): 
    """
    TODO: 
    """
    substrates = [s.strip() for s in substrates_model.split(';')]
    query_lower = substrate_api.lower()
    for substrate in substrates:
        if substrate.lower() == query_lower:
            return substrate
    return None


def fuzzy_match_substrate(substrate_api, substrates_model, threshold=90):
    """
    TODO: 
    """
    substrates = [s.strip() for s in substrates_model.split(';')]
    
    best_match, score, _ = process.extractOne(substrate_api, substrates, scorer=fuzz.token_sort_ratio)
    
    if score >= threshold:
        return best_match, score
    return None, score


# ---------------------------------------------
# Check  
# ---------------------------------------------


def check_enzyme(api_output, kcat_dict):
    """
    TODO: Write the documentation for this function.
    """
    enzyme = kcat_dict['uniprot_model']
    return api_output[api_output['UniProtKB_AC'] == enzyme]


def check_substrate(api_output, kcat_dict, api):
    """
    TODO: Write the documentation for this function.
    """
    if api == 'sabio_rk':
        # 1. Try to match using KEGG reaction ID if provided
        if kcat_dict.get('KEGG_rxn_id'):
            kegg_rxn_matches = api_output[api_output['KeggReactionID'] == kcat_dict['KEGG_rxn_id']]
            if not kegg_rxn_matches.empty:
                result = find_rxn_direction(kcat_dict, kegg_rxn_matches)
                if not result.empty:
                    return result
                
        # 2. If no match or no KEGG, try to match using substrate names
        model_substrates = [s.strip().lower() for s in kcat_dict.get('substrates_name', '').split(';') if s.strip()]
        
        def has_substrate_match(row_substrates):
            api_substrates = [s.strip().lower() for s in row_substrates.split(';') if s.strip()]
            return any(substrate in api_substrates for substrate in model_substrates)
        
        substrate_matches = api_output[api_output['Substrate'].apply(has_substrate_match)]
        return substrate_matches

    elif api == 'brenda':
        # No info on KEGG reaction ID in BRENDA, so we only check substrates
        model_substrates = [s.strip().lower() for s in kcat_dict.get('substrates_name', '').split(';') if s.strip()]
        return api_output[api_output['Substrate'].str.lower().isin(model_substrates)]
    
    else:
        raise ValueError(f"Unknown API: {api}. Supported APIs are 'sabio_rk' and 'brenda'.")


def check_organism(api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    organism = general_criteria["Organism"]
    return api_output[api_output["Organism"] == organism]
    

def check_temp_and_ph(api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    temp_min, temp_max = general_criteria["Temperature"]
    ph_min, ph_max = general_criteria["pH"]

    return api_output[
        (api_output["Temperature"].between(temp_min, temp_max)) &
        (api_output["pH"].between(ph_min, ph_max))
    ]


def check_variant(api_output, general_criteria):
    """
    TODO: Write the documentation for this function.
    """
    variant = general_criteria["Enzyme Variant"]
    return api_output[api_output["Enzyme Variant"] == variant]

# ---------------------------------------------
# Matching functions 
# ---------------------------------------------


def match_exact(kcat_dict, api_output, general_criteria):
    api_output = check_enzyme(api_output, kcat_dict)
    if api_output.empty:
        return None, None
    api_output = check_temp_and_ph(api_output, general_criteria)
    if api_output.empty:
        return None, None
    api_output = check_organism(api_output, general_criteria)
    if api_output.empty:
        return None, None
    api_output = check_substrate(api_output, kcat_dict, api='sabio_rk')
    if api_output.empty:
        return None, None
    return create_kcat_value(api_output), 1



def match_relax_enzyme(kcat_dict, api_output, general_criteria, api):
    api_output = check_temp_and_ph(api_output, general_criteria)
    if api_output.empty:
        return None, None
    api_output = check_organism(api_output, general_criteria)
    if api_output.empty:
        return None, None
    api_output = check_substrate(api_output, kcat_dict, api)
    if api_output.empty:
        return None, None
    return create_kcat_value(api_output), 2


def match_relax_temp_pH(kcat_dict, api_output, general_criteria, api):
    api_output = check_organism(api_output, general_criteria)
    if api_output.empty:
        return None, None
    api_output = check_substrate(api_output, kcat_dict, api)
    if api_output.empty:
        return None, None
    return create_kcat_value(api_output), 3


def match_relax_organism(kcat_dict, api_output, api):
    api_output = check_substrate(api_output, kcat_dict, api)
    if api_output.empty:
        return None, None
    return create_kcat_value(api_output), 4


def match_relax_substrate(api_output):
    return create_kcat_value(api_output), 5


# ---------------------------------------------
# Find best match
# ---------------------------------------------


def find_best_match(kcat_dict, api_output, general_criteria, api):
    """
    TODO: Write the documentation for this function.
    Attempt to find the best match in a hierarchical way:
    1. Perfect match
    2. Relax enzyme
    3. Relax pH and temperature
    4. Relax organism
    5. Relax substrate
    999. Should not happen
    """
    # 1. Perfect match (only for SABIO-RK)
    if api == 'sabio_rk':
        kcat, matching_value = match_exact(kcat_dict, api_output, general_criteria)
        if kcat is not None:
            return kcat, matching_value
    # 2. Relax the enzyme
    kcat, matching_value = match_relax_enzyme(kcat_dict, api_output, general_criteria, api)
    if kcat is not None:
        return kcat, matching_value
    # 3. Relax temperature and pH
    kcat, matching_value = match_relax_temp_pH(kcat_dict, api_output, general_criteria, api)
    if kcat is not None:
        return kcat, matching_value
    # 4. Relax organism
    kcat, matching_value = match_relax_organism(kcat_dict, api_output, api)
    if kcat is not None:
        return kcat, matching_value
    # 5. Relax substrate
    kcat, matching_value = match_relax_substrate(api_output)
    if kcat is not None:
        return kcat, matching_value
    # 6. No match found
    logging.warning('%s: Matching failed' % (kcat_dict['ec_code'])) # Should not happen
    return None, 999


# ---------------------------------------------
# Generate kcat value
# ---------------------------------------------


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
    
#     # For SABIO-RK 
#     api_output = pd.read_csv('in_progress/api_output_test/sabio_rk_test.tsv', sep='\t')
#     kcat, matching_value = find_best_match(kcat_dict, api_output, general_criteria, api='sabio_rk')
#     print(f"Matched kcat: {kcat}, Matching value: {matching_value}")

#     # For BRENDA
#     api_output = pd.read_csv('in_progress/api_output_test/brenda_test.tsv', sep='\t')
#     kcat, matching_value = find_best_match(kcat_dict, api_output, general_criteria, api='brenda')
#     print(f"Matched kcat: {kcat}, Matching value: {matching_value}")