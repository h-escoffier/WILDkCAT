import pandas as pd


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
            "Direction": row["Direction"],
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