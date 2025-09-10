import logging
import pandas as pd
import numpy as np

from ..machine_learning.catapro import create_catapro_input_file, integrate_catapro_predictions
from ..utils.generate_reports import report_prediction_input
from ..utils.manage_warnings import DedupFilter


# TODO: Add a warning if there is the same SMILE for multiple KEGG IDs ?


# --- Format ---


def format_output(kcat_df, limit_matching_score):
    """
    Formats the kcat DataFrame by selecting and rounding kcat values based on matching score and prediction availability.

    Parameters:
        kcat_df (pandas.DataFrame) : Input DataFrame containing kcat-related columns and prediction results.
        limit_matching_score (float) : Threshold for the matching score to determine whether to use predicted kcat values.

    Returns: 
        pandas.DataFrame : Formatted DataFrame with selected and rounded kcat values, reordered columns, and updated source information.
    """
    kcat_df = kcat_df.rename(columns={"kcat": "kcat_source", "kcat_db": "kcat_source_db"})

    def choose_row(row):
        if pd.notna(row["kcat_source"]):
            if row["matching_score"] >= limit_matching_score and pd.notna(row["catapro_predicted_kcat_s"]):
                return pd.Series([row["catapro_predicted_kcat_s"], "catapro"])
            else:
                return pd.Series([row["kcat_source"], row["kcat_source_db"]])
        else:
            if pd.notna(row["catapro_predicted_kcat_s"]):
                return pd.Series([row["catapro_predicted_kcat_s"], "catapro"])
            else:
                return pd.Series([np.nan, np.nan])

    # Add final kcat + db
    kcat_df[["kcat", "kcat_db"]] = kcat_df.apply(choose_row, axis=1)

    # Round numeric columns
    kcat_df["kcat"] = kcat_df["kcat"].round(4)
    kcat_df["kcat_source"] = kcat_df["kcat_source"].round(4)
    kcat_df["kcat_id_percent"] = kcat_df["kcat_id_percent"].round(2)
    kcat_df["catapro_predicted_kcat_s"] = kcat_df["catapro_predicted_kcat_s"].round(1)

    # Reorder columns
    kcat_df = kcat_df[[
        "rxn", "KEGG_rxn_id", "ec_code", "direction", 
        "substrates_name", "substrates_kegg", "products_name", "products_kegg", 
        "genes_model", "uniprot_model",
        "kcat", "kcat_source", "catapro_predicted_kcat_s", 
        "kcat_source_db", "kcat_db", "matching_score",
        "kcat_substrate", "kcat_organism", "kcat_enzyme", "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent"
        ]]

    return kcat_df


# --- Main ---


def run_prediction_part1(kcat_file_path: str, 
                         output_path: str, 
                         limit_matching_score: int, 
                         report: bool = True) -> None:
    """
    Processes kcat data file to generate input files for CataPro prediction.
    Optionally, it can produce a summary report of the processed data.

    Parameters:
        kcat_file_path (str): Path to the input kcat data file.
        output_path (str): Path to save the generated CataPro input CSV file.
        limit_matching_score (int): Threshold for filtering entries based on matching score.
        report (bool, optional): Whether to generate a report using the retrieved data (default: True). 
    """
    # Add a deduplication filter to the logger
    logging.getLogger().addFilter(DedupFilter())

    # Read the kcat file
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')

    # Subset rows with no values or matching score above the limit
    kcat_df = kcat_df[(kcat_df['matching_score'] >= limit_matching_score) | (kcat_df['matching_score'].isnull())]
    # Drop rows with no UniProt ID or no substrates_kegg
    before_duplicates_filter = len(kcat_df) - 1 
    kcat_df = kcat_df[kcat_df['uniprot_model'].notnull() & kcat_df['substrates_kegg'].notnull()]
    nb_missing_enzymes = before_duplicates_filter - len(kcat_df)
    
    # Generate CataPro input file
    catapro_input_df, substrates_to_smiles_df, report_statistics = create_catapro_input_file(kcat_df)

    # Save the CataPro input file and substrates to SMILES mapping
    catapro_input_df.to_csv(output_path, sep=',', index=True)
    substrates_to_smiles_df.to_csv(output_path.replace('.csv', '_substrates_to_smiles.tsv'), sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    # Add statistics 
    report_statistics["missing_enzymes"] = nb_missing_enzymes

    if report:
        report_prediction_input(catapro_input_df, report_statistics)


def run_prediction_part2(kcat_file_path: str, 
                         catapro_predictions_path: str, 
                         substrates_to_smiles_path: str, 
                         output_path: str,
                         limit_matching_score: int) -> None:
    """
    Runs the second part of the kcat prediction pipeline by integrating Catapro predictions,
    mapping substrates to SMILES, formatting the output, and optionally generating a report.
    
    Parameters:
        kcat_file_path (str): Path to the input kcat TSV file.
        catapro_predictions_path (str): Path to the Catapro predictions CSV file.
        substrates_to_smiles_path (str): Path to the TSV file mapping substrates to SMILES.
        output_path (str): Path to save the formatted output TSV file.
        limit_matching_score (float): Threshold for taking predictions over retrieved values.
    """ 
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    substrates_to_smiles = pd.read_csv(substrates_to_smiles_path, sep='\t')
    catapro_predictions_df = pd.read_csv(catapro_predictions_path, sep=',')
    kcat_df = integrate_catapro_predictions(kcat_df, 
                                            substrates_to_smiles,
                                            catapro_predictions_df
                                            )
    
    # Save the output as a TSV file
    kcat_df = format_output(kcat_df, limit_matching_score)
    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")


# if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))
    
    # Test : Main function

    # run_prediction_part1("output/ecoli_kcat_brenda.tsv", -1, "output/machine_learning/ecoli_catapro_input.csv")
    # run_prediction_part2("output/ecoli_kcat_brenda.tsv", 
    #                      "output/machine_learning/ecoli_catapro_output.csv", 
    #                      "output/machine_learning/ecoli_catapro_input_substrates_to_smiles.tsv", 
    #                      8, 
    #                      "output/ecoli_kcat_full.tsv")
    
    # run_prediction_part1("output/yeast_kcat_brenda.tsv", -1, "output/machine_learning/yeast_catapro_input.csv")
    # run_prediction_part2("output/yeast_kcat_brenda.tsv", 
    #                      "output/machine_learning/yeast_catapro_output.csv", 
    #                      "output/machine_learning/yeast_catapro_input_substrates_to_smiles.tsv", 
    #                      8, 
    #                      "output/yeast_kcat_full.tsv")