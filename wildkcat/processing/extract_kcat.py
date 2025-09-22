import re
import os 
import logging
import datetime 
import pandas as pd
from tqdm import tqdm
from functools import lru_cache
from cobra.io import load_json_model, load_matlab_model, read_sbml_model

from ..api.api_utilities import safe_requests_get, retry_api
from ..api.uniprot_api import identify_catalytic_enzyme
from ..utils.manage_warnings import DedupFilter
from ..utils.generate_reports import report_extraction


# --- Load Model ---


def read_model(model_path: str):
    """
    Reads a metabolic model from a given path.
    
    Parameters:
        model_path (str): Path to a model file.

    Returns:
        model (COBRA.Model): The COBRA model object.
    """
    if model_path.endswith(".json"):
        return load_json_model(model_path)
    elif model_path.endswith(".mat"):
        return load_matlab_model(model_path)
    elif model_path.endswith((".xml", ".sbml")):
        return read_sbml_model(model_path)
    else:
        logging.error(f"Unsupported model file format: {model_path}")


# --- KEGG API --- 


@lru_cache(maxsize=None)
def is_ec_code_transferred(ec_code):
    """
    Checks if a given EC code has been transferred according to the KEGG database.

    Parameters:
        ec_code (str): The EC code to check (e.g., '1.1.1.1').

    Returns:
        bool or None: 
            - True if the EC code has been transferred.
            - False if the EC code has not been transferred.
            - None if the KEGG API request fails.

    Logs:
        A warning if the EC code has been transferred.
    """
    url = f'https://rest.kegg.jp/list/{ec_code}'
    safe_get_with_retry = retry_api(max_retries=4, backoff_factor=2)(safe_requests_get)
    response = safe_get_with_retry(url)
    if not response:
        return None
    if "Transferred to" in response.text:
        logging.warning(f"EC code {ec_code} transferred to {response.text.split('Transferred to', 1)[1].lower().strip()}")
        return True
    return False


# --- Generate kcat output ---


def parse_gpr(gpr_str):
    """
    Parses a Gene-Protein-Reaction (GPR) rule string into a list of gene groups.

    The function interprets 'or' as separating alternative gene groups,
    and 'and' as combining genes within a group. Parentheses are used
    to group genes for 'and' relationships.

    Parameters:
        gpr_str (str): The GPR rule string, e.g., '(gene1 and gene2) or gene3'.

    Returns:
        List[List[str]]: A list of gene groups, where each group is a list of gene names.
                         Each inner list represents genes that must all be present ('and'),
                         and each outer list represents alternative gene sets ('or').
    """
    if not gpr_str:
        return []
    or_groups = re.split(r'\s+or\s+', gpr_str, flags=re.IGNORECASE)
    parsed_groups = [
        [g.strip() for g in re.split(r'\s+and\s+', group.replace("(", "").replace(")", ""), flags=re.IGNORECASE) if g.strip()]
        for group in or_groups
    ]
    return [g for g in parsed_groups if g]


def split_metabolites(metabolites):
    """
    Splits a dictionary of metabolites into two lists: metabolite names and their corresponding KEGG IDs.

    Parameters:
        metabolites (dict): A dictionary where keys are metabolite objects and values are their coefficients.
                            Each metabolite object is expected to have 'name', 'id', and 'annotation' attributes.

    Returns:
        tuple: Two lists:
            - names (list of str): Names or IDs of substrate metabolites (where coefficient < 0).
            - kegg_ids (list of str): Corresponding KEGG compound IDs for each substrate metabolite. If not available, an empty string is used.
    """
    names, kegg_ids = [], []
    for m, coeff in metabolites.items():
        if coeff < 0:  # substrate
            name = m.name if m.name else m.id
            kegg = m.annotation.get("kegg.compound")
            if isinstance(kegg, list):
                kegg = kegg[0]
            names.append(name)
            kegg_ids.append(kegg if kegg else "")
    return names, kegg_ids


def create_kcat_output(model):
    """
    Generates a DataFrame summarizing kcat-related information for each reaction in a metabolic model.
    For each reaction, the function extracts KEGG reaction IDs, EC codes, substrates, products, gene-protein-reaction (GPR) associations and UniProt IDs.
    It processes both forward and reverse directions for reversible reactions. 
    The resulting DataFrame contains one row per unique combination of EC code, substrates, products (direction) and gene group.

    Parameters: 
        model (cobra.Model) : A metabolic model. 

    Returns:
        df (pandas.DataFrame) : A DataFrame with columns including reaction ID, KEGG reaction ID, EC code, direction, substrate/product names and KEGG IDs, model genes, UniProt IDs, KEGG genes, and intersection genes.
        report_statistics (dict) : A dictionary with statistics for the report, including the number of incomplete/incorrect EC codes and EC for which kcat values were transferred.
    """
    rows = []
    set_incomplete_ec_codes, set_transferred_ec_codes = set(), set()
    ec_pattern = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
    
    for rxn in tqdm(model.reactions, desc=f"Processing {model.id} reactions"):
        kegg_rxn_id = rxn.annotation.get("kegg.reaction")
        if isinstance(kegg_rxn_id, list):
            kegg_rxn_id = ";".join(kegg_rxn_id)
        ec_codes = rxn.annotation.get("ec-code")
        if not ec_codes:
            continue
        if isinstance(ec_codes, str):
            ec_codes = [ec_codes]
        subs_names, subs_keggs, prod_names, prod_keggs = [], [], [], []
        for m, coeff in rxn.metabolites.items():
            name = m.name if m.name else m.id
            kegg = m.annotation.get("kegg.compound")
            if isinstance(kegg, list):
                kegg = kegg[0]
            if coeff < 0:
                subs_names.append(name)
                subs_keggs.append(kegg if kegg else "")
            elif coeff > 0:
                prod_names.append(name)
                prod_keggs.append(kegg if kegg else "")
        gpr_groups = parse_gpr(rxn.gene_reaction_rule)
        for ec in ec_codes:

            if not ec_pattern.match(ec):
                if ec not in set_incomplete_ec_codes:
                    set_incomplete_ec_codes.add(ec)
                    logging.warning(f"EC code {ec} is not in the correct format")
                continue
            
            is_transferred = is_ec_code_transferred(ec)
            if is_transferred or is_transferred == None:
                if ec not in set_transferred_ec_codes:
                    set_transferred_ec_codes.add(ec)
                    # No continue to be able to identify the number of rxn, kcat dropped due to inconsistency
            
            if not gpr_groups:
                for direction, sn, sk, pn, pk in [
                    ("forward", subs_names, subs_keggs, prod_names, prod_keggs),
                    ("reverse", prod_names, prod_keggs, subs_names, subs_keggs)
                ] if rxn.reversibility else [("forward", subs_names, subs_keggs, prod_names, prod_keggs)]:
                    rows.append({
                        "rxn": rxn.id,
                        "rxn_kegg": kegg_rxn_id,
                        "ec_code": ec,
                        "direction": direction,
                        "substrates_name": ";".join(sn),
                        "substrates_kegg": ";".join(sk),
                        "products_name": ";".join(pn),
                        "products_kegg": ";".join(pk),
                        "genes": "",
                        "uniprot": "",
                        "catalytic_enzyme": "",
                        "warning": ""
                    })
                continue

            for genes_group in gpr_groups:
                genes_group = [g.strip() for g in genes_group if g.strip()]
                uniprot_ids = []
                for gene in genes_group:
                    try:
                        uniprot = model.genes.get_by_id(gene).annotation.get("uniprot")
                        if uniprot:
                            if isinstance(uniprot, list):
                                uniprot_ids.extend(uniprot)
                            else:
                                uniprot_ids.append(uniprot)
                    except KeyError:
                        continue
                uniprot_ids = list(set(uniprot_ids))
                
                # Identify catalytic enzyme
                if len(uniprot_ids) > 1:
                    catalytic_enzyme = identify_catalytic_enzyme(";".join(uniprot_ids), ec)
                else: 
                    catalytic_enzyme = uniprot_ids[0] if uniprot_ids else None

                # Warning 
                warning = (
                    "none" if catalytic_enzyme is None
                    else "" if len(catalytic_enzyme.split(';')) == 1
                    else "multiple"
                )

                for direction, sn, sk, pn, pk in [
                    ("forward", subs_names, subs_keggs, prod_names, prod_keggs),
                    ("reverse", prod_names, prod_keggs, subs_names, subs_keggs)
                ] if rxn.reversibility else [("forward", subs_names, subs_keggs, prod_names, prod_keggs)]:
                    rows.append({
                        "rxn": rxn.id,
                        "rxn_kegg": kegg_rxn_id,
                        "ec_code": ec,
                        "direction": direction,
                        "substrates_name": ";".join(sn),
                        "substrates_kegg": ";".join(sk),
                        "products_name": ";".join(pn),
                        "products_kegg": ";".join(pk),
                        "genes": ";".join(genes_group),
                        "uniprot": ";".join(uniprot_ids),
                        "catalytic_enzyme": catalytic_enzyme, 
                        "warning": warning
                    })

    # Create output 
    df = pd.DataFrame(rows).drop_duplicates(
        subset=["ec_code", "genes", "substrates_name", "products_name", "direction"]
    )

    before_ec_filter = len(df)
    before_nb_reactions = df['rxn'].nunique()
    df = df[~df["ec_code"].isin(set_transferred_ec_codes)]
    nb_lines_dropped = before_ec_filter - len(df)
    nb_rxns_dropped = before_nb_reactions - df['rxn'].nunique()
    
    logging.info("Total of possible kcat values: %d", len(df))

    report_statistics = {
        "incomplete_ec_codes": len(set_incomplete_ec_codes),
        "transferred_ec_codes": len(set_transferred_ec_codes),
        "nb_of_reactions_due_to_unconsistent_ec": nb_rxns_dropped,
        "nb_of_lines_dropped_due_to_unconsistent_ec": nb_lines_dropped
    }

    return df, report_statistics


# --- Main ---


def run_extraction(model_path: str, 
                   output_path: str, 
                   report: bool = True) -> None:
    """
    Extracts kcat-related data from a metabolic model and generates output files and an optional HTML report.

    Parameters:
        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).
        output_path (str): Path to the output file (TSV format).
        report (bool, optional): Whether to generate an HTML report (default: True).
    """
    # Intitialize logging
    os.makedirs("logs", exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"logs/extract_{timestamp}.log"
    logging.getLogger().addFilter(DedupFilter())
    logging.basicConfig(filename=filename, encoding='utf-8', level=logging.INFO)
    
    # Run extraction
    model = read_model(model_path)
    df, report_statistics = create_kcat_output(model)

    df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")
    if report:
        report_extraction(model, df, report_statistics)
