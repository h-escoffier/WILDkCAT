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


# if __name__ == "__main__":
#     print("Starting kcat extraction from model...")
#     extract_kcat_from_model(
#         model_path='data/model/Human-GEM.xml',
#         output_file='output/Human-GEM_kcat.tsv'
#     )
#     print("Kcat extraction completed and saved to output/Human-GEM_kcat.tsv")