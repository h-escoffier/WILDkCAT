import re 
import logging
import requests
import pandas as pd 
from tqdm import tqdm
from functools import lru_cache
from cobra.io import load_json_model, load_matlab_model, read_sbml_model


def read_model(model_path: str):
    """
    Reads a metabolic model from a given path.
    
    Parameters:
    - model_path: str
        Path to a model file.
    
    Returns:
    - model: COBRA.Model
        The COBRA model object.
    """
    # Check the file extension
    if model_path.endswith(".json"):
        return load_json_model(model_path)
    elif model_path.endswith(".mat"):
        return load_matlab_model(model_path)
    elif model_path.endswith(".xml") or model_path.endswith(".sbml"):
        return read_sbml_model(model_path)
    else:
        raise ValueError(f"Unsupported model file format: {model_path}")


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
    ec_url = f"https://rest.kegg.jp/link/{organism_code}/enzyme:{ec_code}"
    response = requests.get(ec_url)
    
    if response.status_code != 200:
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

    for rxn in tqdm(model.reactions, desc=f"Processing {model.name} reactions"):
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
    return df


def run_extraction(model_path, organism_code, report=True):
    """
    Extracts kcat-related data from a metabolic model and generates output files and an optional HTML report.
    Parameters:
        model_path (str): Path to the metabolic model file 
        organism_code (str): Organism code of the model (e.g., 'eco' for E. coli).
        report (bool, optional): If True, generates an HTML report with summary statistics. Defaults to True.

    """
    model = read_model(model_path)
    logging.info(f"Model loaded with {len(model.reactions)} reactions.")
    df = create_kcat_output(model, organism_code)
    df.to_csv("output/ecoli_kcat.tsv", sep='\t', index=False)
    logging.info("Output saved to 'output/ecoli_kcat.tsv'")
    
    if report: 

        # Model statistics
        num_model_reactions = len(model.reactions)
        num_model_metabolites = len(model.metabolites)
        num_model_genes = len(model.genes)

        # Kcat extraction statistics
        num_reactions = df['rxn'].nunique()
        num_ec_codes = df['ec_code'].nunique()
        num_kegg_rxn_ids = df['KEGG_rxn_id'].nunique()
        num_ec_kegg_pairs = df[['ec_code', 'KEGG_rxn_id']].drop_duplicates().shape[0]

        # Coverage statistics
        rxn_coverage = 100.0 * num_reactions / num_model_reactions if num_model_reactions else 0

        html = f"""
        <html>
        <head><title>Kcat Extraction Report</title></head>
        <body>
        <h1>Kcat Extraction Report</h1>
        <h2>Model Overview</h2>
        <ul>
            <li><b>Model name:</b> {model.name}</li>
            <li><b>Number of reactions:</b> {num_model_reactions}</li>
            <li><b>Number of metabolites:</b> {num_model_metabolites}</li>
            <li><b>Number of genes:</b> {num_model_genes}</li>
        </ul>
        <h2>Kcat Extraction Statistics</h2>
        <ul>
            <li><b>Number of reactions with kcat informations:</b> {num_reactions} ({rxn_coverage:.1f}% of model reactions)</li>
            <li><b>Number of unique EC codes:</b> {num_ec_codes}</li>
            <li><b>Number of unique KEGG reaction IDs:</b> {num_kegg_rxn_ids}</li>
            <li><b>Number of unique (EC code, KEGG rxn ID) pairs:</b> {num_ec_kegg_pairs}</li>
            <li><b>Total rows in output:</b> {len(df)}</li>
        </ul>
        </body>
        </html>
        """

        # Save report
        report_path = "reports/kcat_summary.html"
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(html)
        logging.info(f"HTML report saved to '{report_path}'")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("Starting EC code extraction...")
    run_extraction("model/e_coli_core.json", 'eco')
    print("EC code extraction completed.")