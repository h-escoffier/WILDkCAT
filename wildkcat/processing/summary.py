import logging
import pandas as pd
import numpy as np

from ..processing.extract_kcat import read_model
from ..utils.generate_reports import report_final


def generate_summary_report(model_path, kcat_file_path):
    """
    TODO: Write the documentation 
    """
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    model = read_model(model_path)
    report_final(model, kcat_df)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Test : Main function 
    generate_summary_report('model/yeast-GEM.xml', "output/yeast_kcat_full.tsv")