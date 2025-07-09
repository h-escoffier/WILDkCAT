import pandas as pd
import requests
import re
from tqdm import tqdm
from time import sleep

# Cache dictionaries
kegg_to_sid_cache = {}
sid_to_cid_cache = {}
cid_to_smiles_cache = {}
name_to_smiles_cache = {}


def get_pubchem_sid(kegg_id):
    """
    Fetches the PubChem Substance ID (SID) corresponding to a KEGG compound ID.

    Parameters:
        kegg_id (str): KEGG compound ID (e.g., 'C00031').

    Returns:
        str: The PubChem SID if found, otherwise None.
    """
    if kegg_id in kegg_to_sid_cache:
        return kegg_to_sid_cache[kegg_id]

    url = f"https://rest.kegg.jp/get/{kegg_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None

    match = re.search(r'PubChem:\s*(\d+)', response.text)
    sid = match.group(1) if match else None
    kegg_to_sid_cache[kegg_id] = sid
    return sid


def sid_to_cid(sid):
    """
    Converts a PubChem Substance ID (SID) to the corresponding Compound ID (CID).

    Parameters:
        sid (str): PubChem Substance ID.

    Returns:
        int or None: The corresponding PubChem Compound ID (CID), or None if not found.
    """
    if sid in sid_to_cid_cache:
        return sid_to_cid_cache[sid]

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/JSON"
    response = requests.get(url)
    cid = None
    if response.status_code == 200:
        try:
            cid = response.json()['InformationList']['Information'][0]['CID'][0]
        except (KeyError, IndexError):
            cid = None
    sid_to_cid_cache[sid] = cid
    return cid


def cid_to_smiles(cid):    
    """
    Converts a PubChem Compound ID (CID) to its corresponding SMILES representation.

    Parameters:
        cid (str): PubChem Compound ID.

    Returns:
        list or None: A list of SMILES strings if found, otherwise None.
    """
    if cid in cid_to_smiles_cache:
        return cid_to_smiles_cache[cid]

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/smiles/txt"
    try:
        response = requests.get(url)
        response.raise_for_status()
        smiles = response.text.strip().split('\n')
        cid_to_smiles_cache[cid] = smiles
        return smiles
    except:
        return None


def name_to_smiles(name):
    """
    Retrieves the SMILES representation(s) for a compound using its common name via the PubChem API.
    Uses a local cache to avoid redundant API requests.

    Parameters:
        name (str): The common name of the compound (e.g., 'glucose').

    Returns:
        list or None: A list of SMILES strings if found, otherwise None.
    """
    if name in name_to_smiles_cache:
        return name_to_smiles_cache[name]

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/smiles/txt"
    try:
        response = requests.get(url)
        response.raise_for_status()
        smiles = response.text.strip().split('\n')
        name_to_smiles_cache[name] = smiles
        return smiles
    except:
        return None


def get_smiles(substance_name, kegg_id):
    """
    Retrieves the SMILES representation of a chemical compound from PubChem.

    Parameters:
        substance_name (str): The compound identifier (either KEGG ID or common name).
        type (str): Type of identifier ('kegg' or 'name').

    Returns:
        list or None: A list of SMILES strings if found, otherwise None.
    """
    # 1. Try KEGG → SID → CID → SMILES
    if pd.notna(kegg_id) and kegg_id.startswith('C') and kegg_id[1:].isdigit():
        sid = get_pubchem_sid(kegg_id)
        if sid:
            cid = sid_to_cid(sid)
            if cid:
                smiles = cid_to_smiles(cid)
                if smiles:
                    return kegg_id, smiles

    # 2. Fallback to name → SMILES
    if pd.notna(substance_name):
        smiles = name_to_smiles(substance_name)
        if smiles:
            return substance_name, smiles

    return None, None


def add_smiles_to_file(kcat_path, output_path, max_entries=None):
    """
    Adds SMILES annotations to a TSV file containing substrate names and KEGG IDs.

    Parameters:
        kcat_path (str): Path to the input TSV file containing 'Substrates_Name' and 'Substrates_KEGG' columns.
        output_path (str): Path to save the output TSV file with a new 'SMILES' column.
        max_entries (int, optional): Maximum number of unique substrates to process. Defaults to None (no limit).

    Returns:
        pandas.DataFrame: The annotated DataFrame with SMILES added.
    """
    df = pd.read_csv(kcat_path, sep='\t')
    smiles_map = {}
    multiple_smiles = []
    not_found = []

    unique_substrates = df[['Substrates_Name', 'Substrates_KEGG']].drop_duplicates()

    for i, row in tqdm(unique_substrates.iterrows(), total=len(unique_substrates), desc="Fetching SMILES"):
        if max_entries and i >= max_entries:
            break

        name = row['Substrates_Name']
        kegg = row['Substrates_KEGG']

        key, smiles = get_smiles(name, kegg)

        if smiles is None:
            not_found.append(key or f"{name}/{kegg}")
            continue

        if len(smiles) == 1:
            smiles_map[key] = smiles[0]
        else:
            smiles_map[key] = smiles[0]  # Choose first for mapping -> Not ideal
            multiple_smiles.append((key, smiles))

        sleep(0.2)

    def resolve_smiles(row):
        kegg = row['Substrates_KEGG']
        name = row['Substrates_Name']
        return smiles_map.get(kegg) or smiles_map.get(name)

    df['SMILES'] = df.apply(resolve_smiles, axis=1)
    df.to_csv(output_path, sep='\t', index=False)

    print(f"\nFile with SMILES saved to: {output_path}")

    if multiple_smiles:
        df_multi = pd.DataFrame(multiple_smiles, columns=["Substrate", "SMILES_options"])
        df_multi.to_csv(output_path.replace('.tsv', '_multiple_smiles.tsv'), sep='\t', index=False)
        print(f"{len(multiple_smiles)} substrates had multiple SMILES options.")

    if not_found:
        df_nf = pd.DataFrame(not_found, columns=["Substrate"])
        df_nf.to_csv(output_path.replace('.tsv', '_not_found.tsv'), sep='\t', index=False)
        print(f"{len(not_found)} substrates not found in PubChem.")

    return df


# Get UniProt sequences for each enzymes 
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


if __name__ == "__main__":
    # Example usage
    print('start')
    add_smiles_to_file(
        kcat_path="output/CataPro/kcat_Human-GEM.tsv",
        output_path="output/CataPro/kcat_Human-GEM_smiles.tsv", 
        max_entries=10
    )
    print('end')