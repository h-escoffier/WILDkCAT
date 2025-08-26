import logging
import pandas as pd
import numpy as np

from wildkcat.processing.extract_kcat import read_model
from wildkcat.utils.generate_reports import report_final


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Test : Main function 
    # kcat_df = pd.read_csv("output/ecoli_kcat_full.tsv", sep='\t')
    # model = read_model('model/e_coli_core.json')
    kcat_df = pd.read_csv("output/yeast_kcat_full.tsv", sep='\t')
    model = read_model('model/yeast-GEM.xml')
    report_final(model, kcat_df)