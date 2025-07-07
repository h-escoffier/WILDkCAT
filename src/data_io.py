import requests
import pandas as pd
from tqdm import tqdm
from time import sleep
import pandas as pd


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
    for uid in tqdm(iterable=uniprot_ids, desc='Fetching Uniprot sequences'):
        sequence_map[uid] = get_uniprot_sequence(uid)
        sleep(0.01)
    df['Sequence'] = df['Uniprot'].map(sequence_map)
    df.to_csv(output_path, sep='\t', index=False)
    print(f"File with sequences saved to: {output_path}")
    return df


def get_pubchem_cid(kegg_id):
    """
    Fetch the PubChem CID corresponding to a KEGG compound ID.
    
    Parameters:
        kegg_id (str): KEGG compound ID (e.g., 'C00031')
    
    Returns:
        str: PubChem CID, or None if not found
    """
    url = f"https://rest.kegg.jp/get/{kegg_id}"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Error fetching KEGG entry for {kegg_id}")
        return None
    lines = response.text.splitlines()
    in_dblinks = False
    for line in lines:
        if line.startswith("DBLINKS"):
            in_dblinks = True
        elif in_dblinks:
            if line.startswith(" "):
                if "PubChem:" in line:
                    return line.split("PubChem:")[1].strip()
            else:
                break 
    print(f"No PubChem CID found for {kegg_id}")
    return None


def get_smiles(substance_name, type):
    if type == 'kegg':
        cid = get_pubchem_cid(substance_name)
        if cid is None:
            return None
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/smiles/txt"
    elif type == 'name':
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{substance_name}/property/smiles/txt"
    try:
        response = requests.get(url)
        response.raise_for_status()
        smiles = response.text.strip().split('\n')
        return smiles
    except requests.exceptions.HTTPError as _:
        print(f"[404] Substance '{substance_name}' not found in PubChem.")
    except Exception as e:
        print(f"An error occurred while fetching '{substance_name}': {e}")
    return None 

def add_smiles_to_file(kcat_path, output_path, max_entries=None):
    df = pd.read_csv(kcat_path, sep='\t')
    smiles_map = {}
    multiple_smiles = []
    not_found = []

    unique_substrates = df[['Substrates_Name', 'Substrates_KEGG']].drop_duplicates()

    for i, row in tqdm(unique_substrates.iterrows(), total=len(unique_substrates), desc='Fetching SMILES'):
        if max_entries and i >= max_entries:
            break

        name = row['Substrates_Name']
        kegg = row['Substrates_KEGG']
        smiles = None
        key = None

        # 1. Try KEGG if valid
        if pd.notna(kegg) and isinstance(kegg, str) and kegg.startswith('C') and kegg[1:].isdigit():
            smiles = get_smiles(kegg, type='kegg')
            key = kegg

        # 2. Fallback to name if KEGG failed or was NaN
        if smiles is None and pd.notna(name):
            smiles = get_smiles(name, type='name')
            key = name

        if smiles is None:
            not_found.append(key or f"{name}/{kegg}")
            continue

        if len(smiles) == 1:
            smiles_map[key] = smiles[0]
        else:
            smiles_map[key] = smiles[0]  # default to first, log ambiguity
            multiple_smiles.append((key, smiles))

        sleep(0.01)

    # Map SMILES back to original DataFrame
    def resolve_smiles(row):
        kegg = row['Substrates_KEGG']
        name = row['Substrates_Name']
        if pd.notna(kegg) and kegg in smiles_map:
            return smiles_map[kegg]
        elif name in smiles_map:
            return smiles_map[name]
        return None

    df['SMILES'] = df.apply(resolve_smiles, axis=1)
    df.to_csv(output_path, sep='\t', index=False)

    print(f"\nFile with SMILES saved to: {output_path}")

    if multiple_smiles:
        df_multi = pd.DataFrame(multiple_smiles, columns=["Substrate", "SMILES_options"])
        multi_path = output_path.replace(".tsv", "_multiple_smiles.tsv")
        df_multi.to_csv(multi_path, sep='\t', index=False)
        print(f"Multiple SMILES found for {len(multiple_smiles)} substrates (saved to {multi_path})")

    if not_found:
        df_not_found = pd.DataFrame(not_found, columns=["Substrate"])
        nf_path = output_path.replace(".tsv", "_not_found.tsv")
        df_not_found.to_csv(nf_path, sep='\t', index=False)
        print(f"{len(not_found)} substrates not found in PubChem (saved to {nf_path})")

    return df


def format_file_for_catapro(kcat_path, output_file):
    """
    Format an enzyme file by extracting and renaming specific columns.
    
    Parameters:
    - input_file (str): Path to the original TSV file.
    - output_file (str): Path to save the formatted TSV file.
    """
    # Load the original file
    df = pd.read_csv(kcat_path, sep='\t')
    # test with the 10 first lines of the df 
    df = df.head(10)
    # Create a new DataFrame with required columns
    formatted_df = pd.DataFrame({
        'Enzyme_id': df['Uniprot'],
        'type': 'wild',
        'sequence': df['Sequence'],
        'smiles': df['SMILES']
    })
    formatted_df.to_csv(output_file, sep='\t', index=False)


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
    # TurNiP
    # df_with_sequences = add_sequences_to_file(
    #     kcat_path="output/Human-GEM_kcat_clean.tsv", 
    #     output_path="output/Human-GEM_kcat_seq.tsv"
    # )
    # remove_nan_values(
    #     kcat_path='output/Human-GEM_kcat_seq.tsv',
    #     output_path='output/Human-GEM_kcat_seq_clean.tsv'
    # )
    # CataPro 
    # add_sequences_to_file(
    #     kcat_path="output/CataPro/kcat_Human-GEM.tsv", 
    #     output_path="output/CataPro/kcat_Human-GEM_seq.tsv"
    # )
    # add_smiles_to_file(kcat_path="output/CataPro/kcat_Human-GEM_seq.tsv",
    #                    output_path="output/CataPro/kcat_Human-GEM_seq_smiles.tsv", 
    #                    max_entries=100)
    format_file_for_catapro(
        kcat_path="output/CataPro/kcat_Human-GEM_seq_smiles.tsv", 
        output_file="output/CataPro/kcat_Human-GEM_formatted.tsv"
    ) 
    print("end")
