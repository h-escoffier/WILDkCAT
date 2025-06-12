from cobra.io import read_sbml_model
import pandas as pd 


def extract_kcat_from_model(model_path, output_file):
    """
    Extracts kinetic parameters from a GEM model and saves them to a TSV file.
    """
    model = read_sbml_model(model_path)
    records = []
    for rxn in model.reactions:
        ec_code = rxn.annotation.get("ec-code", None)
        if ec_code is None:
            continue

        ec_list = ec_code if isinstance(ec_code, list) else [ec_code]

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

        # KEGG compound separation
        substrates, products = [], []
        for met, coeff in rxn.metabolites.items():
            kegg_id = met.annotation.get("kegg.compound", "NaN")
            if isinstance(kegg_id, list):
                kegg_id = kegg_id[0]
            if coeff < 0:
                substrates.append(kegg_id)
            else:
                products.append(kegg_id)

        substrate_str = ";".join(substrates) if substrates else "NaN"
        product_str = ";".join(products) if products else "NaN"

        for ec in ec_list:
            for up in uniprot_ids:
                records.append({
                    "RxnID": rxn.id,
                    "EC_Code": ec,
                    "Uniprot": up,
                    "KEGG_Substrate": substrate_str,
                    "KEGG_Product": product_str
                })

    # Write to TSV
    df = pd.DataFrame(records)
    df.to_csv(output_file, sep="\t", index=False)


def remove_nan_values(kcat_path, output_path=None):
    """
    Removes rows where any column contains the string 'NaN'.
    """
    import pandas as pd

    df = pd.read_csv(kcat_path, sep="\t", dtype=str)  # Force strings
    # Keep rows where none of the values are the string "NaN"
    df_cleaned = df[~df.apply(lambda row: row.astype(str).str.contains("NaN")).any(axis=1)]

    if output_path:
        df_cleaned.to_csv(output_path, sep="\t", index=False)
    else:
        df_cleaned.to_csv(kcat_path, sep="\t", index=False)


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


# if __name__ == "__main__":
#     print("start")
#     extract_kcat_from_model(
#         model_path='data/model/Human-GEM.xml',
#         output_file='output/Human-GEM_kcat.tsv'
#     )
#     remove_nan_values(
#         kcat_path='output/Human-GEM_kcat.tsv',
#         output_path='output/Human-GEM_kcat_cleaned.tsv'
#     )
#     generate_enzyme_xlsx_with_uniprot(
#         tsv_path='output/Human-GEM_kcat_cleaned.tsv',
#         xlsx_output_path='output/Human-GEM_kcat.xlsx',
#         max_rows=500
#     )
#     print("end")
