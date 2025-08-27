import pandas as pd 
import requests
from time import sleep
from Bio import Entrez
from tqdm import tqdm


def taxon_id_to_name(taxon_id, email="hugues.escoffier@uni.lu"):
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="taxonomy", id=str(taxon_id), retmode="xml")
        records = Entrez.read(handle)
        return records[0]["ScientificName"]
    except Exception as e:
        print(f"Error retrieving taxon name for ID {taxon_id}: {e}")
        return None


def get_taxon_id(uniprot_id):
    """Query UniProt API to get the taxon ID of a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data.get("organism", {}).get("taxonId")
        else:
            print(f"Failed to retrieve {uniprot_id}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return None


def filter_by_taxon_uniprot_api(df: pd.DataFrame, taxon_id: int, sleep_time: float = 0.01) -> pd.DataFrame:
    """Filter a DataFrame by matching Uniprot ID to the given taxon ID using UniProt API."""

    uniprot_ids = df['Uniprot ID'].unique().tolist()

    # Map Uniprot ID to taxon ID
    id_to_taxon = {}
    for uid in tqdm(iterable=uniprot_ids):
        tax_id = get_taxon_id(uid)
        if tax_id is None:
            print(f"Taxon ID not found for {uid}.")
        id_to_taxon[uid] = tax_id
        sleep(sleep_time) 

    # Add Taxon ID column
    df = df.copy()
    df['Taxon ID'] = df['Uniprot ID'].map(id_to_taxon)

    # Filter by desired taxon ID
    return df[df['Taxon ID'] == taxon_id]  


def retrieve_kcat_from_organism(brenda_kcat, sabio_kcat, uniprot_kcat, organism_taxon_id):
    scientific_name = str(taxon_id_to_name(organism_taxon_id))
    # Extract kcat values for a specific organism in BRENDA  
    brenda_df = pd.read_csv(brenda_kcat, sep='\t')
    brenda_df_filtered = brenda_df[brenda_df['ORGANISM'] == scientific_name]
    print(brenda_df_filtered.shape)

    # Extract kcat values for a specific organism in SABIO-RK
    sabio_df = pd.read_csv(sabio_kcat, sep='\t')
    sabio_df_filtered = filter_by_taxon_uniprot_api(sabio_df, organism_taxon_id)
    print(sabio_df_filtered.shape)

    # Extract kcat values for a specific organism in UniProt
    uniprot_df = pd.read_csv(uniprot_kcat, sep='\t')
    uniprot_df_filtered = filter_by_taxon_uniprot_api(uniprot_df, organism_taxon_id)
    print(uniprot_df_filtered.shape)
    print(uniprot_df_filtered.head())

    # combined dataframe based on Kcat values, Uniprot ID, substrate and product


retrieve_kcat_from_organism(
    brenda_kcat='data/kcat/BRENDA_kcat.tsv',
    sabio_kcat='data/kcat/Sabio_kcat.tsv',
    uniprot_kcat='data/kcat/Uniprot_kcat.tsv',
    organism_taxon_id=9606)

