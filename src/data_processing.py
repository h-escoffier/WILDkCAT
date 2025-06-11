import pickle
import pandas as pd

def load_pckl_data(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj

kcat_data = load_pckl_data('data/kcat/final_kcat_dataset.pkl')
kcat_data.to_csv('data/kcat/final_kcat_dataset.tsv', sep='\t', index=False)