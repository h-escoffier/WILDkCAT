import os 
import logging
import pandas as pd
import numpy as np

from ..processing.extract_kcat import read_model
from ..utils.generate_reports import report_final


def generate_summary_report(model_path: str,
                            output_folder: str) -> None:
    """
    Generate a HTML report summarizing the kcat extraction, retrieval and prediction for a given model. 

    Parameters:
        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).
        output_folder (str): Path to the output folder where the kcat file is located.
    """
    # Read the kcat file
    if not os.path.exists(output_folder):
        raise FileNotFoundError(f"The specified output folder '{output_folder}' does not exist.")
    
    kcat_full_file_path = os.path.join(output_folder, "kcat_full.tsv")
    kcat_retrieve_file_path = os.path.join(output_folder, "kcat_retrieved.tsv")
    if os.path.isfile(kcat_full_file_path):
        kcat_df = pd.read_csv(kcat_full_file_path, sep='\t')
        model = read_model(model_path)
        report_final(model, kcat_df, output_folder)
    elif os.path.isfile(kcat_retrieve_file_path):
        logging.warning(f"The file 'kcat_full.tsv' is not present in the folder '{output_folder}' the general report will be done without predicted values.")
        model = read_model(model_path)
        kcat_df = pd.read_csv(kcat_retrieve_file_path, sep='\t')
        report_final(model, kcat_df, output_folder)
    else: 
        raise FileNotFoundError(f"The specified folder '{output_folder}' does not contain the files: 'kcat_full.tsv', 'kcat_retrieve.tsv'. Please run at least the extraction step.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Test : Main function 
    generate_summary_report('model/yeast-GEM.xml', "output/yeast_kcat_full.tsv")