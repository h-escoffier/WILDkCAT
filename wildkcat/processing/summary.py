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
    
    kcat_file_path = os.path.join(output_folder, "kcat_full.tsv")
    if not os.path.isfile(kcat_file_path):
        raise FileNotFoundError(f"The specified file '{kcat_file_path}' does not exist in the output folder. Please run the full pipeline.")
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    model = read_model(model_path)
    report_final(model, kcat_df, output_folder)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Test : Main function 
    generate_summary_report('model/yeast-GEM.xml', "output/yeast_kcat_full.tsv")