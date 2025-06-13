from cobra.io import read_sbml_model
import pandas as pd 


def extract_kcat_from_model(model_path, output_file):
    """
    Extracts kinetic parameters from a GEM model and saves them to a TSV file,
    including both KEGG and ChEBI metabolite annotations.
    
    Parameters:
    - model_path: Path to the SBML model file.
    - output_file: Path to the output TSV file.
    """
    model = read_sbml_model(model_path)
    records = []

    for rxn in model.reactions:
        ec_code = rxn.annotation.get("ec-code", None)
        if ec_code is None:
            continue

        ec_list = ec_code if isinstance(ec_code, list) else [ec_code]

        # Extract UniProt IDs from associated genes
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

        # Extract metabolite annotations
        kegg_substrates, kegg_products = [], []
        chebi_substrates, chebi_products = [], []

        for met, coeff in rxn.metabolites.items():
            kegg = met.annotation.get("kegg.compound", "NaN")
            if isinstance(kegg, list):
                kegg = kegg[0]
            elif not isinstance(kegg, str):
                kegg = str(kegg)

            chebi = met.annotation.get("chebi", "NaN")
            if isinstance(chebi, list):
                chebi = chebi[0]
            elif not isinstance(chebi, str):
                chebi = str(chebi)

            if coeff < 0:
                kegg_substrates.append(kegg)
                chebi_substrates.append(chebi)
            else:
                kegg_products.append(kegg)
                chebi_products.append(chebi)

        for ec in ec_list:
            for up in uniprot_ids:
                records.append({
                    "RxnID": rxn.id,
                    "EC_Code": ec,
                    "Uniprot": up,
                    "kegg.compound_Substrate": ";".join(kegg_substrates) if kegg_substrates else "NaN",
                    "kegg.compound_Product": ";".join(kegg_products) if kegg_products else "NaN",
                    "chebi_Substrate": ";".join(chebi_substrates) if chebi_substrates else "NaN",
                    "chebi_Product": ";".join(chebi_products) if chebi_products else "NaN"
                })

    df = pd.DataFrame(records)
    df.to_csv(output_file, sep="\t", index=False)


def remove_nan_values(kcat_path, output_path=None):
    """
    Filters out rows where:
    - Both KEGG and ChEBI metabolite annotations contain 'NaN'
    - Uniprot field contains 'NaN'

    Keeps only rows where either KEGG or ChEBI metabolite annotations are fully valid.
    If both are valid, prefers KEGG. Outputs columns:
    RxnID, EC_Code, Uniprot, Substrates, Products
    """

    df = pd.read_csv(kcat_path, sep="\t", dtype=str).fillna("NaN")

    def is_valid(field):
        """Check if all values in a semicolon-separated field are valid (i.e., not 'NaN')."""
        return all(x.strip() != "NaN" for x in str(field).split(";"))

    cleaned_records = []

    for _, row in df.iterrows():
        if row["Uniprot"] == "NaN":
            continue  # Skip rows without a valid Uniprot ID

        kegg_sub_valid = is_valid(row["kegg.compound_Substrate"])
        kegg_prod_valid = is_valid(row["kegg.compound_Product"])
        chebi_sub_valid = is_valid(row["chebi_Substrate"])
        chebi_prod_valid = is_valid(row["chebi_Product"])

        if kegg_sub_valid and kegg_prod_valid:
            substrates = row["kegg.compound_Substrate"]
            products = row["kegg.compound_Product"]
        elif chebi_sub_valid and chebi_prod_valid:
            substrates = row["chebi_Substrate"]
            products = row["chebi_Product"]
        else:
            continue  # Skip if neither set is valid

        cleaned_records.append({
            "RxnID": row["RxnID"],
            "EC_Code": row["EC_Code"],
            "Uniprot": row["Uniprot"],
            "Substrates": substrates,
            "Products": products
        })

    cleaned_df = pd.DataFrame(cleaned_records)

    if output_path:
        cleaned_df.to_csv(output_path, sep="\t", index=False)
    else:
        cleaned_df.to_csv(kcat_path, sep="\t", index=False)


def generate_enzyme_xlsx_with_uniprot(tsv_path, xlsx_output_path, max_rows=None):
    """
    Generate an XLSX file with UniProt ID in 'Enzyme', and KEGG compound IDs in 'Substrates' and 'Products'.

    Parameters:
    - tsv_path (str): Input TSV file with columns: RxnID, EC Code, Uniprot, KEGG substrate, KEGG product
    - xlsx_output_path (str): Output XLSX file path
    - max_rows (int, optional): Number of rows to limit processing (for testing)
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    
    rows = []
    for i, row in df.iterrows():
        if max_rows and i >= max_rows:
            break
        rows.append({
            "Enzyme": row["Uniprot"],
            "Substrates": row["KEGG_Substrate"],
            "Products": row["KEGG_Product"]
        })

    out_df = pd.DataFrame(rows)
    out_df.to_excel(xlsx_output_path, index=False)
    print(f"XLSX file saved to: {xlsx_output_path}")


if __name__ == "__main__":
    print("start")
    extract_kcat_from_model(
        model_path='model/Human-GEM.xml',
        output_file='output/Human-GEM_kcat.tsv'
    )
    remove_nan_values(
        kcat_path='output/Human-GEM_kcat.tsv',
        output_path='output/Human-GEM_kcat_clean.tsv'
    )
    # generate_enzyme_xlsx_with_uniprot(
    #     tsv_path='output/Human-GEM_kcat_clean.tsv',
    #     xlsx_output_path='output/Human-GEM_kcat.xlsx',
    #     max_rows=500
    # )
    print("end")
