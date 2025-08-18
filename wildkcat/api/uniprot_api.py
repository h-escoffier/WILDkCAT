import requests
from functools import lru_cache


# --- UniProt API ---


@lru_cache(maxsize=None)
def convert_uniprot_to_sequence(uniprot_id):
    """
    Convert a UniProt accession ID to its corresponding amino acid sequence.

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
        # logging.warning(f"Failed to retrieve sequence for UniProt ID {uniprot_id}")
        return None


if __name__ == "__main__":
    # Test : Send a request to UniProt API
    uniprot_id = "Q16774"
    seq = convert_uniprot_to_sequence(uniprot_id)
    print(seq)