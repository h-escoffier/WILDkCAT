import pickle
import pandas as pd 


# conda create -n 'env_pkl' python=3.9 pandas=1.3 suds-jurko
# conda activate env_pkl

# Pickle files were download from: https://zenodo.org/records/8367052


def load_pckl_data(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


# kcat_data = load_pckl_data('../data/kcat/Uniprot_kcat.pkl')
# kcat_data.to_csv('../data/kcat/Uniprot_kcat.tsv', sep='\t', index=False)
# kcat_data = load_pckl_data('../data/kcat/BRENDA_kcat.pkl')
# kcat_data.to_csv('../data/kcat/BRENDA_kcat.tsv', sep='\t', index=False)
# kcat_data = load_pckl_data('../data/kcat/Sabio_kcat.pkl')
# kcat_data.to_csv('../data/kcat/Sabio_kcat.tsv', sep='\t', index=False)


def format_sabio_kcat_data(sabio_kcat = 'data/kcat/Sabio_kcat.tsv'): 
    sabio_df = pd.read_csv(sabio_kcat, sep='\t')
    sabio_df = sabio_df.copy()
    # Split the 'Uniprot ID' column by space into lists
    sabio_df['Uniprot ID'] = sabio_df['Uniprot ID'].str.split()
    # Explode the DataFrame so each UniProt ID gets its own row
    sabio_df = sabio_df.explode('Uniprot ID').reset_index(drop=True)
    sabio_df.to_csv('Sabio_kcat.tsv', sep='\t', index=False)
