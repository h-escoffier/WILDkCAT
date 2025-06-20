import requests
import pandas as pd
from tqdm import tqdm
from time import sleep


def get_uniprot_sequence(uniprot_id):
    """
    Fetches the amino acid sequence for a given UniProt ID.

    Parameters:
        uniprot_id (str): The UniProt accession ID.

    Returns:
        str: The amino acid sequence, or None if not found.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        fasta = response.text
        lines = fasta.splitlines()
        sequence = ''.join(lines[1:])  # Skip the header
        return sequence
    else:
        print(f"Failed to fetch sequence for UniProt ID {uniprot_id}. HTTP status code: {response.status_code}")
        return None


def add_sequences_to_file(kcat_path, output_path):
    """
    Adds UniProt sequences to the kcat DataFrame and saves the result to a file.

    Parameters:
        kcat_path (str): Path to the input kcat file.
        output_path (str): Path to save the output file with sequences.

    Returns:
        pd.DataFrame: DataFrame containing the kcat data with added sequences.
    """
    df = pd.read_csv(kcat_path, sep='\t')

    uniprot_ids = df['Uniprot'].unique()
    sequence_map = {}

    for uid in tqdm(iterable=uniprot_ids):
        sequence_map[uid] = get_uniprot_sequence(uid)
        sleep(0.01)
    df['Sequence'] = df['Uniprot'].map(sequence_map)

    df.to_csv(output_path, sep='\t', index=False)
    print(f"File with sequences saved to: {output_path}")

    return df


def remove_nan_values(kcat_path, output_path=None):
    """
    Removes rows with NaN values from the kcat DataFrame.

    Parameters:
        kcat_path (str): Path to the input kcat file.
        output_path (str, optional): Path to save the output file without NaN values.

    Returns:
        pd.DataFrame: DataFrame without NaN values.
    """
    df = pd.read_csv(kcat_path, sep='\t')
    df_cleaned = df.dropna()
    if output_path:
        df_cleaned.to_csv(output_path, sep='\t', index=False)
    return df_cleaned


if __name__ == "__main__":
    print("start")
    df_with_sequences = add_sequences_to_file(
        kcat_path="output/Human-GEM_kcat_clean.tsv", 
        output_path="output/Human-GEM_kcat_seq.tsv"
    )
    remove_nan_values(
        kcat_path='output/Human-GEM_kcat_seq.tsv',
        output_path='output/Human-GEM_kcat_seq_clean.tsv'
    )
    print("end")
