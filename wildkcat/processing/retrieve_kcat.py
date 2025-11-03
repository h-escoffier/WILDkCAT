import os
import datetime
import time
import logging
from dotenv import load_dotenv
import numpy as np
import pandas as pd
from tqdm import tqdm
from functools import lru_cache 
import pandas as pd

from ..api.sabio_rk_api import get_turnover_number_sabio, get_enzyme_sabio
from ..api.brenda_api import get_turnover_number_brenda, get_enzyme_brenda

from ..utils.matching import find_best_match
from ..utils.manage_warnings import DedupFilter
from ..utils.generate_reports import report_retrieval


@lru_cache(maxsize=None)
def get_turnover_number(ec_code, database='both'): 
    """
    Retrieves turnover number (kcat) data from specified enzyme databases and returns a merged DataFrame.

    Parameters: 
        enzyme_uniprot (str): The UniProt ID of the enzyme.
        database (str, optional): Specifies which database(s) to query for kcat values. 
            Options are 'both' (default), 'brenda', or 'sabio_rk'.

    Returns: 
        pd.DataFrame: A DataFrame containing kcat data from the selected database(s), with columns unified across sources.

    Raises:
        ValueError: If an invalid database option is provided.
    """
    df_brenda = pd.DataFrame()
    df_sabio = pd.DataFrame()

    if database in ('both', 'brenda'): 
        df_brenda = get_turnover_number_brenda(ec_code)
    if database in ('both', 'sabio_rk'):
        df_sabio = get_turnover_number_sabio(ec_code)
        time.sleep(1)  
    if database not in ('both', 'brenda', 'sabio_rk'):
        raise ValueError("Invalid database option. Choose from 'both', 'brenda', or 'sabio_rk'.")

    # Get columns 
    all_columns = set(df_brenda.columns).union(df_sabio.columns)

    # Merge all outputs
    df_brenda = df_brenda.reindex(columns=all_columns, fill_value=None)
    df_sabio = df_sabio.reindex(columns=all_columns, fill_value=None)
    non_empty_dfs = [df for df in [df_brenda, df_sabio] if not df.empty]
    if non_empty_dfs:
        df = pd.concat(non_empty_dfs, ignore_index=True)
    else:
        df = pd.DataFrame(columns=list(all_columns))
    return df


@lru_cache(maxsize=None)
def get_enzyme(enzyme_uniprot, organism, database='both'):
    """
    Retrieves enzyme data from specified databases and returns a merged DataFrame.

    Parameters: 
        enzyme_uniprot (str): The UniProt ID of the enzyme.
        database (str, optional): Specifies which database(s) to query for enzyme data. 
            Options are 'both' (default), 'brenda', or 'sabio_rk'.

    Returns: 
        pd.DataFrame: A DataFrame containing kcat data from the selected database(s), with columns unified across sources.

    Raises:
        ValueError: If an invalid database option is provided.
    """
    df_brenda = pd.DataFrame()
    df_sabio = pd.DataFrame()

    if database in ('both', 'brenda'): 
        df_brenda = get_enzyme_brenda(enzyme_uniprot, organism)
    if database in ('both', 'sabio_rk'):
        df_sabio = get_enzyme_sabio(enzyme_uniprot)
        time.sleep(1)  
    
    # Get columns 
    all_columns = set(df_brenda.columns).union(df_sabio.columns)

    # Merge all outputs
    df_brenda = df_brenda.reindex(columns=all_columns, fill_value=None)
    df_sabio = df_sabio.reindex(columns=all_columns, fill_value=None)
    non_empty_dfs = [df for df in [df_brenda, df_sabio] if not df.empty]
    if non_empty_dfs:
        df = pd.concat(non_empty_dfs, ignore_index=True)
    else:
        df = pd.DataFrame(columns=list(all_columns))
    return df


def extract_kcat(kcat_dict, general_criteria, database='both'): 
    """
    Extracts the best matching kcat value from a given set of criteria.

    Parameters:
        kcat_dict (dict): Dictionary containing enzyme information.
        general_criteria (dict): Dictionary specifying matching criteria.
        database (str, optional): Specifies which database(s) to query for kcat values. 
            Options are 'both' (default), 'brenda', or 'sabio_rk'.

    Returns:
        tuple: 
            - best_candidate (dict or None): The best matching kcat entry, or None if no match is found.
            - best_score (int or float): The score of the best candidate, or 17 if no match is found in the database.
    """
    if kcat_dict['ec_code'] == '' or kcat_dict['warning_ec'] in ['incomplete', 'transferred']:
        api_output = get_enzyme(kcat_dict['catalytic_enzyme'], general_criteria['Organism'], database) 
    else: 
        api_output = get_turnover_number(kcat_dict['ec_code'], database)
    
    if api_output.empty: 
        return None, 17
            
    best_score, best_candidate = find_best_match(kcat_dict, api_output, general_criteria)
    return best_candidate, best_score


def merge_ec(kcat_df: pd.DataFrame):
    """
    Merge entries with the same combination of reaction and substrate that differ only in EC numbers.
    Select the kcat entry with the highest matching score; in case of a tie, use the following priorities:
    1. Highest matching_score
    2. Highest sequence_score
    3. Closest organism_score
    4. Highest kcat_value

    Parameters:
        kcat_df (pd.DataFrame): DataFrame containing kcat data.  

    Returns:
        update_kcat_df (pd.DataFrame): Updated DataFrame with merged EC numbers.
    """
    
    # Sort by the selection criteria
    kcat_df_sorted = kcat_df.sort_values(
        by=['matching_score', 'kcat_id_percent', 'kcat_organism_score', 'kcat'],
        ascending=[True, False, True, False]
    )

    # Handle missing ec
    kcat_df['ec_code'] = kcat_df['ec_code'].fillna('').astype(str)

    def merge_ec_codes(x):
        ec_list = [ec for ec in x if ec and ec.lower() != 'nan']
        return ';'.join(sorted(set(ec_list))) if ec_list else 'NA'
    
    # Merge EC numbers for each reaction-substrate pair
    ec_merged = (
        kcat_df
        .groupby(['rxn', 'substrates_name', 'products_kegg', 'genes', 'uniprot'])['ec_code']
        .apply(merge_ec_codes)
        .rename('ec_codes')
    )


    best_entries = kcat_df_sorted.groupby(['rxn', 'substrates_name', 'products_kegg', 'genes', 'uniprot'], as_index=False).first()

    # Add merged EC numbers to best entries
    update_kcat_df = best_entries.merge(ec_merged, on=['rxn', 'substrates_name', 'products_kegg', 'genes', 'uniprot'])

    # Reorder columns to place 'ec_codes' next to 'ec_code'
    update_kcat_df = update_kcat_df[
        [
            'rxn', 'rxn_kegg', 'ec_code', 'ec_codes', 'direction',
            'substrates_name', 'substrates_kegg', 'products_name', 'products_kegg',
            'genes', 'uniprot', 'catalytic_enzyme', 'warning_ec', 'warning_enz',
            'kcat', 'matching_score', 'kcat_substrate', 'kcat_organism', 'kcat_enzyme',
            'kcat_temperature', 'kcat_ph', 'kcat_variant', 'kcat_db',
            'kcat_id_percent', 'kcat_organism_score'
        ]
    ]

    return update_kcat_df


def run_retrieval(output_folder: str,
                  organism: str,
                  temperature_range: tuple,
                  pH_range: tuple,
                  database: str = 'both',
                  report: bool = True) -> None:
    """
    Retrieves closest kcat values from specified databases for entries in a kcat file, applies filtering criteria, 
    and saves the results to an output file.
    
    Parameters:
        output_folder (str): Path to the output folder where the results will be saved.
        organism (str): Organism scientific name (e.g. "Escherichia coli", "Homo sapiens").
        temperature_range (tuple): Acceptable temperature range for filtering (min, max).
        pH_range (tuple): Acceptable pH range for filtering (min, max).
        database (str, optional): Specifies which database(s) to query for kcat values. 
            Options are 'both' (default), 'brenda', or 'sabio_rk'.
        report (bool, optional): Whether to generate an HTML report using the retrieved data (default: True).        
    """
    # Load environment variables
    load_dotenv()

    # Create a dict with the general criterias
    general_criteria = {
        "Organism": organism,
        "Temperature": temperature_range,
        "pH": pH_range
    }

    # Read the kcat file
    if not os.path.exists(output_folder):
        raise FileNotFoundError(f"The specified output folder '{output_folder}' does not exist.")
    
    # Intitialize logging
    os.makedirs(os.path.join(output_folder, "logs"), exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"logs/retrieval_{timestamp}.log"
    logging.getLogger().addFilter(DedupFilter())
    logging.basicConfig(filename=os.path.join(output_folder, filename), encoding='utf-8', level=logging.INFO)
    
    kcat_file_path = os.path.join(output_folder, "kcat.tsv")
    if not os.path.isfile(kcat_file_path):
        raise FileNotFoundError(f"The specified file '{kcat_file_path}' does not exist in the output folder. Please run the function 'run_extraction()' first.")
    
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    
    # Initialize new columns
    kcat_df['kcat'] = None
    kcat_df['matching_score'] = None

    # Add data of the retrieve kcat values
    kcat_df['kcat_substrate'] = None
    kcat_df['kcat_organism'] = None
    kcat_df['kcat_enzyme'] = None
    kcat_df['kcat_temperature'] = None
    kcat_df['kcat_ph'] = None
    kcat_df['kcat_variant'] = None
    kcat_df['kcat_db'] = None

    # Retrieve kcat values from databases
    request_count = 0
    for row in tqdm(kcat_df.itertuples(), total=len(kcat_df), desc="Retrieving kcat values"):
        kcat_dict = row._asdict()
        
        # Extract kcat and matching score
        best_match, matching_score = extract_kcat(kcat_dict, general_criteria, database=database)
        kcat_df.loc[row.Index, 'matching_score'] = matching_score

        request_count += 1
        if request_count % 300 == 0:
            time.sleep(10)
        
        if best_match is not None:
            # Assign results to the main dataframe
            kcat_df.loc[row.Index, 'kcat'] = best_match['adj_kcat']
            kcat_df.loc[row.Index, 'kcat_substrate'] = best_match['Substrate']
            kcat_df.loc[row.Index, 'kcat_organism'] = best_match['Organism']
            kcat_df.loc[row.Index, 'kcat_enzyme'] = best_match['UniProtKB_AC']
            kcat_df.loc[row.Index, 'kcat_temperature'] = best_match['adj_temp']
            kcat_df.loc[row.Index, 'kcat_ph'] = best_match['pH']
            kcat_df.loc[row.Index, 'kcat_variant'] = best_match['EnzymeVariant']
            kcat_df.loc[row.Index, 'kcat_db'] = best_match['db']
            kcat_df.loc[row.Index, 'kcat_id_percent'] = best_match['id_perc']
            kcat_df.loc[row.Index, 'kcat_organism_score'] = best_match['organism_score']

    # Select only one kcat value per reaction and substrate
    kcat_df = merge_ec(kcat_df)

    output_path = os.path.join(output_folder, "kcat_retrieved.tsv")
    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    if report:
        report_retrieval(kcat_df, output_folder)


if __name__ == "__main__":
    # # Test : Send a request for a specific EC number
    # kcat_dict = {
    #     'ec_code': '1.16.3.1',
    #     'rxn_kegg': 'R00078',
    #     'uniprot': 'P36649',
    #     'catalytic_enzyme': 'P36649',
    #     'substrates_name': 'Fe2+ mitochondria;H+;O2 O2', 
    #     'warning_ec': ''
    # }

    # general_criteria ={
    #     'Organism': 'Escherichia coli', 
    #     'Temperature': (20, 40), 
    #     'pH': (6.5, 7.5)
    # }

    # output = extract_kcat(kcat_dict, general_criteria, database='both')
    # print(output)

    # Test : Run the retrieve function

    # run_retrieval(
    #     kcat_file_path="output/ecoli_kcat.tsv",
    #     output_path="output/ecoli_kcat_both.tsv",
    #     # output_path="output/ecoli_kcat_sabio.tsv",
    #     organism="Escherichia coli",
    #     temperature_range=(20, 40),
    #     pH_range=(6.5, 7.5),
    #     database='both', 
    #     # database='brenda', 
    #     # database='sabio_rk', 
    #     report=False
    # ) 

    # run_retrieval(
    #     kcat_file_path="output/yeast_kcat.tsv",
    #     output_path="output/yeast_kcat_brenda.tsv",
    #     # output_path="output/yeast_kcat_sabio.tsv",
    #     organism="Saccharomyces cerevisiae",
    #     temperature_range=(18, 38),
    #     pH_range=(4.0, 8.0),
    #     database='brenda', 
    #     # database='sabio_rk', 
    #     report=True
    # ) 

    # Test : Generate report
    # df = pd.read_csv("output/yeast_kcat_brenda.tsv", sep='\t')
    df = pd.read_csv("in_progress/ecoli_v3/kcat_retrieved.tsv", sep='\t')
    report_retrieval(df, output_folder="in_progress/ecoli_v3")