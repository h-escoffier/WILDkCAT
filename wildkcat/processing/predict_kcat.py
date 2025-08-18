import logging
import pandas as pd
from tqdm import tqdm

from wildkcat.machine_learning.catapro import create_catapro_input_file, integrate_catapro_predictions

# TODO: Generates reports


# --- Main ---


def run_catapro_part1(kcat_file_path, limit_matching_score, output_path, report=True):
    """
    TODO 
    """
    # Read the kcat file
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')

    # Subset rows with no values or matching score above the limit
    kcat_df = kcat_df[(kcat_df['matching_score'] > limit_matching_score) | (kcat_df['matching_score'].isnull())]
    # Drop rows with no UniProt ID or no substrates_kegg
    kcat_df = kcat_df[kcat_df['uniprot_model'].notnull() & kcat_df['substrates_kegg'].notnull()]
    
    # Generate CataPro input file
    catapro_input_df, substrates_to_smiles_df, report_statistics = create_catapro_input_file(kcat_df)

    logging.info(
    f"Generated CataPro input file with {len(catapro_input_df)} entries.\n"
    f"Statistics summary:\n"
    f"  - Reactions covered: {report_statistics['rxn_covered']}\n"
    f"  - Cofactors identified: {report_statistics['cofactor_identified']}\n"
    f"  - Reactions with multiple UniProt IDs: {report_statistics['multiple_uniprot']}\n"
    f"  - KEGG entries without match: {report_statistics['kegg_no_matching']}"
)

    # Save the CataPro input file and substrates to SMILES mapping
    catapro_input_df.to_csv(output_path, sep=',', index=True)
    substrates_to_smiles_df.to_csv(output_path.replace('.csv', '_substrates_to_smiles.tsv'), sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    if report:
        pass 


def run_catapro_part2(kcat_file_path, catapro_predictions_path, substrates_to_smiles_path, output_path, report=True):
    """
    TODO: Write the documentation
    """ 
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    substrates_to_smiles = pd.read_csv(substrates_to_smiles_path, sep='\t')
    catapro_predictions_df = pd.read_csv(catapro_predictions_path, sep=',')
    kcat_df = integrate_catapro_predictions(kcat_df, 
                                            substrates_to_smiles,
                                            catapro_predictions_df
                                            )
    
    # Save the output as a TSV file
    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    if report:
        pass # TODO: Implement report generation



if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))

    # Test : Integrate CataPro predictions into kcat file
    # kcat_df = pd.read_csv("output/ecoli_kcat_sabio.tsv", sep='\t')
    # substrates_to_smiles = pd.read_csv('in_progress/ml_test/substrates_to_smiles.tsv', sep='\t')
    # integrate_catapro_predictions(kcat_df, substrates_to_smiles, "in_progress/ml_test/catapro_output.csv", "in_progress/ml_test/ecoli_kcat_catapro.tsv")

    # Test : Main function
    logging.basicConfig(level=logging.INFO)
    # run_catapro_part1("output/ecoli_kcat_test_brenda.tsv", 5, "output/machine_learning/ecoli_catapro_input.csv")
    # run_catapro_part1("output/yeast_kcat_test_brenda.tsv", -1, "output/machine_learning/yeast_catapro_input.csv")
    run_catapro_part2("output/ecoli_kcat_brenda.tsv", 
                      "in_progress/ml_test/catapro_output.csv", 
                      "output/machine_learning/ecoli_catapro_input_substrates_to_smiles.tsv", 
                      "output/ecoli_kcat_catapro.tsv")