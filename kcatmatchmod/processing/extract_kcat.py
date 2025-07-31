import re 
import logging
import requests
import pandas as pd 
from tqdm import tqdm
from functools import lru_cache

from kcatmatchmod.reports.generate_reports import report_extraction
from kcatmatchmod.utils.model_function import read_model


# TODO: Implement Retry mechanism for API calls (to deal with status 502, 503, 504)
# TODO: The KEGG genes have to match the genes in the model to be able to check the intersection
# TODO: The report can be improved with more statistics, e.g. number of kcats dropped because the EC have been transferred etc.


def get_kegg_genes_by_ec(organism_code: str, ec_code: str):
    """
    Retrieve genes for a given organism and EC number from KEGG.

    Parameters:
    - organism_code : str
        KEGG organism code (e.g., 'hsa' for human, 'eco' for E. coli).
    - ec_code : str
        EC number (e.g., '1.1.1.1').

    Returns:
    - list
        List of gene identifiers (KEGG format).
    """
    # Get KEGG orthologs (or reactions) associated with the EC code
    url = f"https://rest.kegg.jp/link/{organism_code}/enzyme:{ec_code}"
    response = requests.get(url)
    
    if response.status_code != 200:
        logging.warning(f"Failed to retrieve data from KEGG API: {response.status_code} - {ec_code}")
        return []  # or raise Exception(f"Failed to retrieve data from KEGG API: {response.status_code}")
    
    if response.text == "\n":
        return []

    genes = []

    for line in response.text.strip().split("\n"):
        _, gene_entry = line.split("\t")
        if gene_entry.startswith(f"{organism_code}:"):
            gene = gene_entry[len(organism_code) + 1:]
            genes.append(gene)

    return genes


@lru_cache(maxsize=None)
def cached_get_kegg_genes_by_ec(organism_code, ec_code):
    """Cached wrapper to avoid repeated KEGG API calls."""
    return get_kegg_genes_by_ec(organism_code, ec_code)


def is_ec_code_transferred(ec_code): 
    """
    TODO: 
    """
    url = f'https://rest.kegg.jp/list/{ec_code}'
    response = requests.get(url)

    if response.status_code != 200:
        logging.error(f"Failed to retrieve data from KEGG API: {response.status_code}")
        return None

    if "Transferred to" in response.text:
        logging.warning(f"EC code {ec_code} has been transferred.")
        return True
    return False


@lru_cache(maxsize=None)
def cached_is_ec_code_transferred(ec_code):
    """Cached wrapper to avoid repeated KEGG API calls."""
    return is_ec_code_transferred(ec_code)


def parse_gpr(gpr_str):
    """
    Parse GPR into groups (OR → multiple entries, AND → combined).
    Example: '(gene1 and gene2) or gene3' ->
             [['gene1','gene2'], ['gene3']]
    """
    if not gpr_str:
        return []

    # Split by 'or' (outer level)
    or_groups = re.split(r'\s+or\s+', gpr_str, flags=re.IGNORECASE)

    parsed_groups = []
    for group in or_groups:
        # Remove parentheses and split 'and'
        genes = re.split(r'\s+and\s+', group.replace("(", "").replace(")", ""), flags=re.IGNORECASE)
        genes = [g.strip() for g in genes if g.strip()]
        if genes:
            parsed_groups.append(genes)

    return parsed_groups


def split_metabolites(metabolites):
    """
    Return two lists: names and KEGG IDs separately.
    """
    names = []
    kegg_ids = []
    for m, coeff in metabolites.items():
        if coeff < 0:  # substrate
            name = m.name if m.name else m.id
            kegg = m.annotation.get("kegg.compound")
            if isinstance(kegg, list):
                kegg = kegg[0]
            names.append(name)
            kegg_ids.append(kegg if kegg else "")
    return names, kegg_ids


def create_kcat_output(model, organism_code):
    rows = []
    transferred_count = 0
    for rxn in tqdm(model.reactions, desc=f"Processing {model.id} reactions"):
        # --- KEGG Reaction ID ---
        kegg_rxn_id = rxn.annotation.get("kegg.reaction")
        if isinstance(kegg_rxn_id, list):
            kegg_rxn_id = ";".join(kegg_rxn_id)

        # --- EC Codes ---
        ec_codes = rxn.annotation.get("ec-code")
        if not ec_codes:
            continue
        if isinstance(ec_codes, str):
            ec_codes = [ec_codes]

        # --- Substrates / Products ---
        subs_names = []
        subs_keggs = []
        prod_names = []
        prod_keggs = []
        for m, coeff in rxn.metabolites.items():
            name = m.name if m.name else m.id
            kegg = m.annotation.get("kegg.compound")
            if isinstance(kegg, list):
                kegg = kegg[0]
            if coeff < 0:  # substrate
                subs_names.append(name)
                subs_keggs.append(kegg if kegg else "")
            elif coeff > 0:  # product
                prod_names.append(name)
                prod_keggs.append(kegg if kegg else "")

        # --- GPR parsing ---
        gpr_groups = parse_gpr(rxn.gene_reaction_rule)

        # --- For each EC code ---
        for ec in ec_codes:
            # KEGG genes for EC
            kegg_genes = cached_get_kegg_genes_by_ec(organism_code, ec)

            if kegg_genes == []:
                is_transferred = cached_is_ec_code_transferred(ec)
                if is_transferred: 
                    transferred_count += 1
                continue

            # If no genes are associated with the EC code in the model 
            if not gpr_groups:
                rows.append({
                    "rxn": rxn.id,
                    "KEGG_rxn_id": kegg_rxn_id,
                    "ec_code": ec,
                    "direction": "forward",
                    "substrates_name": ";".join(subs_names),
                    "substrates_kegg": ";".join(subs_keggs),
                    "products_name": ";".join(prod_names),
                    "products_kegg": ";".join(prod_keggs),
                    "genes_model": "",
                    "kegg_genes": ";".join(kegg_genes),
                    "intersection_genes": ""
                })
                if rxn.reversibility:
                    rows.append({
                        "rxn": rxn.id,
                        "KEGG_rxn_id": kegg_rxn_id,
                        "ec_code": ec,
                        "direction": "reverse",
                        "substrates_name": ";".join(prod_names),
                        "substrates_kegg": ";".join(prod_keggs),
                        "products_name": ";".join(subs_names),
                        "products_kegg": ";".join(subs_keggs),
                        "genes_model": "",
                        "kegg_genes": ";".join(kegg_genes),
                        "intersection_genes": ""
                    })
                continue

            # Otherwise process per GPR group
            for genes_group in gpr_groups:
                # Retieve UniprotIDs for the genes in the group
                genes_group = [g.strip() for g in genes_group if g.strip()]
                uniprot_ids = []
                for gene in genes_group: 
                    uniprot = model.genes.get_by_id(gene).annotation.get("uniprot")
                    if uniprot:
                        if isinstance(uniprot, list):
                            uniprot_ids.extend(uniprot)
                        else:
                            uniprot_ids.append(uniprot)
                uniprot_ids = list(set(uniprot_ids))
                intersection = list(set(genes_group) & set(kegg_genes))

                # Forward row
                rows.append({
                    "rxn": rxn.id,
                    "KEGG_rxn_id": kegg_rxn_id,
                    "ec_code": ec,
                    "direction": "forward",
                    "substrates_name": ";".join(subs_names),
                    "substrates_kegg": ";".join(subs_keggs),
                    "products_name": ";".join(prod_names),
                    "products_kegg": ";".join(prod_keggs),
                    "genes_model": ";".join(genes_group),
                    "uniprot_model": ";".join(uniprot_ids),
                    "kegg_genes": ";".join(kegg_genes),
                    "intersection_genes": ";".join(intersection) if intersection else ""
                })

                # Reverse row if reversible
                if rxn.reversibility:
                    rows.append({
                        "rxn": rxn.id,
                        "KEGG_rxn_id": kegg_rxn_id,
                        "ec_code": ec,
                        "direction": "reverse",
                        "substrates_name": ";".join(prod_names),
                        "substrates_kegg": ";".join(prod_keggs),
                        "products_name": ";".join(subs_names),
                        "products_kegg": ";".join(subs_keggs),
                        "genes_model": ";".join(genes_group),
                        "uniprot_model": ";".join(uniprot_ids),
                        "kegg_genes": ";".join(kegg_genes),
                        "intersection_genes": ";".join(intersection) if intersection else ""
                    })

    # Remove duplicates
    df = pd.DataFrame(rows).drop_duplicates(
        subset=["ec_code", "genes_model", "substrates_name", "products_name", "direction"]
    )

    logging.info("Total of possible kcat values: %d", len(df))
    return df, transferred_count


def run_extraction(model_path, output_path, organism_code, report=True):
    """
    Extracts kcat-related data from a metabolic model and generates output files and an optional HTML report.
    Parameters:
        model_path (str): Path to the metabolic model file 
        organism_code (str): Organism code of the model (e.g., 'eco' for E. coli).
        report (bool, optional): If True, generates an HTML report with summary statistics. Defaults to True.

    """
    model = read_model(model_path)
    logging.info(f"Model loaded with {len(model.reactions)} reactions.")
    df, transferred = create_kcat_output(model, organism_code)
    df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    if report:
        report_extraction(model, df, transferred)
        

if __name__ == "__main__":
    # Test : Run extraction on a model
    logging.basicConfig(level=logging.INFO)
    run_extraction("model/e_coli_core.json", "output/ecoli_kcat.tsv", 'eco')
    # run_extraction("model/Human-GEM.xml", "output/human_gem_kcat.tsv", 'hsa')