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
    For each reaction, the function extracts KEGG reaction IDs, EC codes (if available), substrates, products, 
    gene-protein-reaction (GPR) associations, and UniProt IDs. 
    It processes both forward and reverse directions for reversible reactions. 
    The resulting DataFrame contains one row per unique combination of EC code, substrates, products (direction), and gene group.

    Parameters: 
        model (cobra.Model): A metabolic model.

    Returns:
        df (pandas.DataFrame): DataFrame with columns:
            - rxn, rxn_kegg, ec_code, direction, substrates/products names & KEGG IDs, genes, uniprot, catalytic_enzyme, warning_ec, warning_enz
        report_statistics (dict): A dictionary with statistics for the report, including the number of incomplete/incorrect EC codes and EC for which kcat values were transferred.
    """
    rows = []
    ec_pattern = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

    for rxn in tqdm(model.reactions, desc=f"Processing {model.id} reactions"):
        kegg_rxn_id = rxn.annotation.get("kegg.reaction")
        if isinstance(kegg_rxn_id, list):
            kegg_rxn_id = ";".join(kegg_rxn_id)

        ec_codes = rxn.annotation.get("ec-code")
        if isinstance(ec_codes, str):
            ec_codes = [ec_codes]
        elif not ec_codes:
            ec_codes = [""]

        # Extract substrates and products
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

        # Parse GPR
        gpr_groups = parse_gpr(rxn.gene_reaction_rule)

        for ec in ec_codes:
            warning_ec = ""

            # Check EC code validity (if present)
            if ec:
                if not ec_pattern.match(ec):
                    warning_ec = "incomplete"
                    logging.warning(f"EC code {ec} is not in the correct format")
                    
                if warning_ec != "incomplete":
                    is_transferred = is_ec_code_transferred(ec)
                    if is_transferred or is_transferred is None:
                        warning_ec = "transferred"

            # If no GPR 
            if not gpr_groups:
                for direction, sn, sk, pn, pk in (
                    [("forward", subs_names, subs_keggs, prod_names, prod_keggs),
                     ("reverse", prod_names, prod_keggs, subs_names, subs_keggs)]
                    if rxn.reversibility else
                    [("forward", subs_names, subs_keggs, prod_names, prod_keggs)]
                ):
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
                        "warning_ec": warning_ec,
                        "warning_enz": "no_gpr"
                    })
                continue

            # Check GPR 
            for genes_group in gpr_groups:
                warning_enz = ""

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
                    catalytic_enzyme = uniprot_ids[0] if uniprot_ids else ""

                # Enzyme-related warnings
                if not catalytic_enzyme:
                    warning_enz = "none"
                elif ";" in catalytic_enzyme:
                    warning_enz = "multiple"

                for direction, sn, sk, pn, pk in (
                    [("forward", subs_names, subs_keggs, prod_names, prod_keggs),
                     ("reverse", prod_names, prod_keggs, subs_names, subs_keggs)]
                    if rxn.reversibility else
                    [("forward", subs_names, subs_keggs, prod_names, prod_keggs)]
                ):
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
                        "warning_ec": warning_ec,
                        "warning_enz": warning_enz
                    })

    # Build final df
    df = pd.DataFrame(rows)

    report_statistics = {
        "nb_missing_ec": (df["ec_code"] == "").sum(),
        "nb_incomplete_ec": (df["warning_ec"] == "incomplete").sum(),
        "nb_transferred_ec": (df["warning_ec"] == "transferred").sum(),
        "nb_missing_gpr": (df["warning_enz"] == "no_gpr").sum(),
        "nb_missing_catalytic_enzyme": (df["warning_enz"] == "none").sum(),
        "nb_multiple_catalytic_enzymes": (df["warning_enz"] == "multiple").sum()
    }

    # Filtering
    rows_before = len(df) - 1 
    rxn_before = df['rxn'].nunique() - 1 

    df = df[~((df["ec_code"] == "") & (df["catalytic_enzyme"] == ""))] # Remove rows without EC and without catalytic enzyme
    df = df[~((df["ec_code"] == "") & (df["warning_enz"] == "none"))]

    rows_exchange = len(df) - 1
    rxn_exchange = df['rxn'].nunique() - 1
    
    # df = df[~df["warning_ec"].isin(["incomplete", "transferred"])] # Remove rows with incomplete or transferred EC codes
    # df.drop(columns=["warning_ec"], inplace=True) # Remove the warning_ec column 

    # rxn_after = df['rxn'].nunique() - 1

    report_statistics.update({
        "nb_of_lines_dropped_no_ec_no_enzyme": rows_before - rows_exchange,
        "nb_of_reactions_dropped_no_ec_no_enzyme": rxn_before - rxn_exchange,
        # "nb_of_reactions_dropped_incomplete_transferred_ec": rxn_exchange - rxn_after
    })

    logging.info("Total of possible kcat values: %d", len(df) - 1)

    return df, report_statistics


# --- Main ---


def run_extraction(model_path: str, 
                   output_folder: str, 
                   report: bool = True) -> None:
    """
    Extracts kcat-related data from a metabolic model and generates output files and an optional HTML report.

    Parameters:
        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).
        output_folder (str): Path to the output folder where all the results will be saved.
        report (bool, optional): Whether to generate an HTML report (default: True).
    """
    # Initialize output folder
    os.makedirs(output_folder, exist_ok=True)

    # Intitialize logging
    os.makedirs(os.path.join(output_folder, "logs"), exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"logs/extract_{timestamp}.log"
    logging.getLogger().addFilter(DedupFilter())
    logging.basicConfig(filename=os.path.join(output_folder, filename), encoding='utf-8', level=logging.INFO)
    
    # Run extraction
    model = read_model(model_path)
    df, report_statistics = create_kcat_output(model)

    # Save output
    output_path = os.path.join(output_folder, "kcat.tsv")
    df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    print(report_statistics)

    if report:
        report_extraction(model, df, report_statistics, output_folder)
