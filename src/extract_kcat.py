import logging
import pandas as pd 
import requests
from tqdm import tqdm

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
    

def retrieve_ec_code(model): 
    """
    Retrieve EC codes from a COBRA model.

    Parameters:
    - model: COBRA.Model
        The COBRA model object.

    Returns:
    - list(ec_codes): list
        A list of unique EC codes found in the model's reactions.
    """
    ec_codes = set()
    for reaction in model.reactions:
        ec = reaction.annotation.get("ec-code")
        if ec:
            ec_codes.update(ec if isinstance(ec, list) else [ec])
    return list(ec_codes)


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

 
def create_kcat_output(model, organism_code):
    """
    Create a DataFrame with detailed reaction info including KEGG reaction IDs,
    EC codes, genes, substrates, and products. Handles reversible reactions by
    including reverse direction as separate entries.

    Parameters:
    - model: COBRA.Model
    - organism_code: str (e.g., 'eco', 'hsa')

    Returns:
    - pd.DataFrame
    """
    rows = []
    for rxn in model.reactions:
        # 1) Reaction ID
        rxn_id = rxn.id
        
        # 2) KEGG reaction annotation (if exists)
        kegg_rxn_id = None
        if "kegg.reaction" in rxn.annotation:
            kegg = rxn.annotation["kegg.reaction"]
            if isinstance(kegg, list):
                kegg_rxn_id = ";".join(kegg)
            else:
                kegg_rxn_id = kegg
        
        # 3) EC codes (may be list or string)
        ec_codes = rxn.annotation.get("ec-code")
        if not ec_codes:
            ec_codes = []
        elif isinstance(ec_codes, str):
            ec_codes = [ec_codes]
        
        # 4) For each EC code, get genes from KEGG
        ec_genes_map = {}
        for ec in ec_codes:
            genes = get_kegg_genes_by_ec(organism_code, ec)
            ec_genes_map[ec] = genes
        
        # 5) Collect genes related to the reaction itself (from model.gene_reaction_rule)
        # This is an alternative gene source, but for now we'll keep focus on KEGG genes.
        # We could also parse rxn.genes or rxn.gene_reaction_rule if needed.
        
        # 6) Prepare substrates and products strings: "met_id:coeff"
        substrates = []
        products = []
        for met, coeff in rxn.metabolites.items():
            # coeff < 0 means substrate, coeff > 0 means product
            s = f"{met.id}({abs(coeff)})"
            if coeff < 0:
                substrates.append(s)
            elif coeff > 0:
                products.append(s)
        
        # 7) Add forward reaction entry
        rows.append({
            "rxn": rxn_id,
            "KEGG_rxn_id": kegg_rxn_id,
            "ec-code": ";".join(ec_codes) if ec_codes else None,
            "genes": ";".join([g for genes in ec_genes_map.values() for g in genes]) if ec_genes_map else None,
            "substrates": ";".join(substrates),
            "products": ";".join(products),
            "direction": "forward"
        })
        
        # 8) If reversible, add reverse reaction entry with substrates and products swapped
        if rxn.reversibility:
            rows.append({
                "rxn": rxn_id,
                "KEGG_rxn_id": kegg_rxn_id,
                "ec-code": ";".join(ec_codes) if ec_codes else None,
                "genes": ";".join([g for genes in ec_genes_map.values() for g in genes]) if ec_genes_map else None,
                "substrates": ";".join(products),   # swapped here
                "products": ";".join(substrates),   # swapped here
                "direction": "reverse"
            })
    
    df = pd.DataFrame(rows)
    return df




def run(model_path, organism_code):
    # model = read_model(model_path)
    # ec_codes = retrieve_ec_code(model)
    # logging.info(f"Retrieved {len(ec_codes)} EC codes from the model.")
    # data = []
    # for ec_code in ec_codes:
    #     genes = get_kegg_genes_by_ec(organism_code, ec_code)
    #     data.append({
    #         "ec_code": ec_code,
    #         "genes": ", ".join(genes) if genes else None
    #     })   
    # # Export results to a CSV file
    # pd.DataFrame(data).to_csv("output/ec_codes.csv", index=False)
    # logging.info("EC codes saved to ec_codes.csv")
    model = read_model(model_path)
    logging.info(f"Model loaded with {len(model.reactions)} reactions.")
    df = create_kcat_output(model, organism_code)
    df.to_csv("output/kcat_reaction_info.csv", index=False)
    logging.info("kcat reaction info saved to output/kcat_reaction_info.csv")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("Starting EC code extraction...")
    run("model/e_coli_core.json", 'eco')
    print("EC code extraction completed.")