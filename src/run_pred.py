from kcat_pred import *
import pandas as pd
import esm


def load_model():
    """
    Load the ESM model and alphabet.
    """
    model, alphabet = torch.hub.load("facebookresearch/esm:v0.4.0", "esm1b_t33_650M_UR50S")
    return model, alphabet


def run_pred(tsv_path, save_path=None): 
    df = pd.read_csv(tsv_path, sep='\t')

    enzymes, substrates, products = [], [], []

    for _, row in df.iterrows():
        enzymes.append(row['Sequence'])
        substrates.append(row['Substrates'])
        products.append(row['Products'])

    df = kcat_prediction(substrates = substrates,
               products = products,
               enzymes = enzymes)

    if save_path:
        df.to_csv(save_path, sep='\t', index=False)
        print(f"Results saved to {save_path}")
    return  df


if __name__ == "__main__":
    print("start")
    load_model()
    run_pred("output/Human-GEM_kcat_clean.tsv", save_path="output/Human-GEM_kcat_pred.tsv")
    print("end")