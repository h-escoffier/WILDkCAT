import logging
import pandas as pd 
from tqdm import tqdm

from cobra.io import load_json_model, load_matlab_model, read_sbml_model


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def read_model(model_path):
    """
    Reads a metabolic model from a given path.
    
    Parameters:
    - model_path: Path to a model file.
    
    Returns:
    - model: The COBRA model object.
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


def is_reversible(rxn, rev_bound=1000):
    """
    Checks if a reaction is reversible.
    
    Parameters:
    - rxn: The reaction object to check.
    - rev_bound: The standard reversibility bound (default is 1000).
    
    Returns:
    - bool: True if the reaction is reversible, False otherwise.
    """
    return rxn.lower_bound == -rev_bound and rxn.upper_bound == rev_bound


def extract_uniprot_ids(gene):
    """Returns a list of UniProt IDs for a gene."""
    ids = gene.annotation.get("uniprot", [])
    if isinstance(ids, str):
        return [ids]
    return list(ids) if ids else []


# Extract values from a GEM model to use TurNiP method 
def extract_kcat_from_model(model_path, output_file):
    """
    Extracts kinetic parameters from a GEM model and saves them to a TSV file,
    including EC codes, UniProt IDs, directionality, and KEGG/ChEBI annotations.
    Reversible reactions are included twice with a 'Direction' column.

    Parameters:
    - model_path: Path to the SBML model file.
    - output_file: Path to the output TSV file.
    """
    model = read_model(model_path)
    records = []

    for rxn in model.reactions:
        ec_code = rxn.annotation.get("ec-code", None)
        if ec_code is None:
            continue
        ec_list = ec_code if isinstance(ec_code, list) else [ec_code]

        # Extract UniProt-annotated genes only
        uniprot_ids = set()
        for gene in rxn.genes:
            if "uniprot" in gene.annotation:
                ids = gene.annotation["uniprot"]
                if isinstance(ids, str):
                    uniprot_ids.add(ids)
                else:
                    uniprot_ids.update(ids)

        if not uniprot_ids:
            uniprot_ids = {"NaN"}

        # Helper to extract metabolite annotations
        def extract_met_annotations():
            kegg_substrates, kegg_products = [], []
            chebi_substrates, chebi_products = [], []

            for met, coeff in rxn.metabolites.items():
                kegg = met.annotation.get("kegg.compound", "NaN")
                kegg = kegg[0] if isinstance(kegg, list) else str(kegg)
                chebi = met.annotation.get("chebi", "NaN")
                chebi = chebi[0] if isinstance(chebi, list) else str(chebi)

                if coeff < 0:
                    kegg_substrates.append(kegg)
                    chebi_substrates.append(chebi)
                else:
                    kegg_products.append(kegg)
                    chebi_products.append(chebi)

            return kegg_substrates, kegg_products, chebi_substrates, chebi_products

        # First direction (forward)
        kegg_sub, kegg_prod, chebi_sub, chebi_prod = extract_met_annotations()

        direction_fwd = "reversible_forward" if rxn.lower_bound == -1000 and rxn.upper_bound == 1000 else "irreversible"

        for ec in ec_list:
            for up in uniprot_ids:
                records.append({
                    "RxnID": rxn.id,
                    "EC_Code": ec,
                    "Uniprot": up,
                    "Direction": direction_fwd,
                    "kegg.compound_Substrate": ";".join(kegg_sub) if kegg_sub else "NaN",
                    "kegg.compound_Product": ";".join(kegg_prod) if kegg_prod else "NaN",
                    "chebi_Substrate": ";".join(chebi_sub) if chebi_sub else "NaN",
                    "chebi_Product": ";".join(chebi_prod) if chebi_prod else "NaN"
                })

        # If reversible rxn add reverse direction 
        if rxn.lower_bound == -1000 and rxn.upper_bound == 1000:
            for ec in ec_list:
                for up in uniprot_ids:
                    records.append({
                        "RxnID": rxn.id,
                        "EC_Code": ec,
                        "Uniprot": up,
                        "Direction": "reversible_reverse",
                        "kegg.compound_Substrate": ";".join(kegg_prod) if kegg_prod else "NaN",
                        "kegg.compound_Product": ";".join(kegg_sub) if kegg_sub else "NaN",
                        "chebi_Substrate": ";".join(chebi_prod) if chebi_prod else "NaN",
                        "chebi_Product": ";".join(chebi_sub) if chebi_sub else "NaN"
                    })

    df = pd.DataFrame(records)
    df.to_csv(output_file, sep="\t", index=False)


def extract_uniprot_from_model(model_path, output_file=None):
    """
    Extracts UniProt IDs and annotated substrates for each reaction in the model.
    Reversible reactions are handled with forward and reverse directions.

    Parameters:
    - model_path (str): Path to the model file.
    - output_file (str, optional): Path to save the output TSV.

    Returns:
    - DataFrame with extracted data.
    """
    model = read_model(model_path)
    all_rxns = []

    for rxn in tqdm(model.reactions, desc=f"Extracting from {model.id}"):
        direction = "reversible" if is_reversible(rxn) else "irreversible"

        # Split substrates based on reaction direction
        substrates_fwd, names_fwd, kegg_fwd = [], [], []
        substrates_rev, names_rev, kegg_rev = [], [], []

        for met, coeff in rxn.metabolites.items():
            kegg = met.annotation.get("kegg.compound", "NaN")
            kegg = kegg[0] if isinstance(kegg, list) else str(kegg)
            if coeff < 0:
                substrates_fwd.append(met.id)
                names_fwd.append(met.name)
                kegg_fwd.append(kegg)
            elif coeff > 0 and direction == "reversible":
                substrates_rev.append(met.id)
                names_rev.append(met.name)
                kegg_rev.append(kegg)

        # Extract UniProt IDs
        uniprot_ids = []
        for gene in rxn.genes:
            uniprot_ids.extend(extract_uniprot_ids(gene))

        for uid in uniprot_ids:
            if direction == "irreversible":
                for idx, sid in enumerate(substrates_fwd):
                    all_rxns.append({
                        "RxnID": rxn.id,
                        "Uniprot": uid,
                        "Direction": direction,
                        "Substrates_ID": sid,
                        "Substrates_Name": names_fwd[idx],
                        "Substrates_KEGG": kegg_fwd[idx]
                    })
            else:
                for idx, sid in enumerate(substrates_fwd):
                    all_rxns.append({
                        "RxnID": rxn.id,
                        "Uniprot": uid,
                        "Direction": "reversible_forward",
                        "Substrates_ID": sid,
                        "Substrates_Name": names_fwd[idx],
                        "Substrates_KEGG": kegg_fwd[idx]
                    })
                for idx, sid in enumerate(substrates_rev):
                    all_rxns.append({
                        "RxnID": rxn.id,
                        "Uniprot": uid,
                        "Direction": "reversible_reverse",
                        "Substrates_ID": sid,
                        "Substrates_Name": names_rev[idx],
                        "Substrates_KEGG": kegg_rev[idx]
                    })

    df = pd.DataFrame(all_rxns)
    if output_file:
        df.to_csv(output_file, sep="\t", index=False)
        logging.info(f"Saved UniProt extraction to {output_file}")
    return df
 

def extract_for_sabio_rk_api(model, output_file=None):
    model = read_model(model)
    results = []

    for rxn in tqdm(model.reactions, desc=f"Extracting from {model.id}"):
        # Retrieve annotations
        kegg_reaction_id = rxn.annotation.get('kegg.reaction') or rxn.annotation.get('kegg_reaction')
        ec_code = rxn.annotation.get('ec-code') or rxn.annotation.get('ec_number')
        
        # Skip if missing
        if not kegg_reaction_id or not ec_code:
            continue

        # Normalize values (list → first element)
        if isinstance(kegg_reaction_id, list):
            kegg_reaction_id = kegg_reaction_id[0]
        if isinstance(ec_code, list):
            ec_code = ec_code[0]

        # Extract substrates and products with KEGG IDs if available
        substrates = []
        products = []
        for met, coeff in rxn.metabolites.items():
            # Try KEGG compound annotation
            kegg_compound = met.annotation.get('kegg.compound') or met.annotation.get('kegg') or None
            if isinstance(kegg_compound, list):
                kegg_compound = kegg_compound[0]

            # If no KEGG, fallback to metabolite ID
            identifier = kegg_compound if kegg_compound else met.id

            # Classify substrate/product
            if coeff < 0:
                substrates.append(identifier)
            elif coeff > 0:
                products.append(identifier)

        results.append({
            'rxn_id': rxn.id,
            'kegg_reaction_id': kegg_reaction_id,
            'ec_code': ec_code,
            'substrates': ";".join(substrates),
            'products': ";".join(products)
        })

    # Convert to DataFrame
    df = pd.DataFrame(results)

    if output_file:
        df.to_csv(output_file, sep="\t", index=False)
        logging.info(f"Saved SABIO-RK extraction to {output_file}")

    return df

    

if __name__ == "__main__":
    print("start")
    # TurNiP
    # extract_kcat_from_model(
    #     model_path='model/Human-GEM.xml',
    #     output_file='output/Human-GEM_kcat.tsv'
    # )
    # remove_nan_values(
    #     kcat_path='output/Human-GEM_kcat.tsv',
    #     output_path='output/Human-GEM_kcat_clean.tsv'
    # )
    # generate_enzyme_xlsx_with_uniprot(
    #     tsv_path='output/Human-GEM_kcat_clean.tsv',
    #     xlsx_output_path='output/Human-GEM_kcat.xlsx',
    #     max_rows=500
    # )
    # CataPro
    # extract_uniprot_from_model(
    #     model_path='model/Human-GEM.xml',
    #     output_file='output/CataPro/kcat_Human-GEM.tsv'
    # ) 
    # extract_uniprot_from_model(
    #     model_path='model/e_coli_core.json',
    #     output_file='output/EColiCore/test.tsv'
    # ) 
    # extract_kcat_from_model(
    #     model_path='model/e_coli_core.json',
    #     output_file='output/EColiCore/test_kcat.tsv'
    # )
    extract_for_sabio_rk_api(
        model='model/Human-GEM.xml', 
        output_file='output/HumanGEM/Human-GEM_sabio_rk.tsv'
    )
    extract_for_sabio_rk_api(
        model='model/e_coli_core.json',
        output_file='output/EColiCore/sabio_rk.tsv'
    )
    print("end")
