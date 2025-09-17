import requests
import logging
from functools import lru_cache


# --- UniProt API ---


@lru_cache(maxsize=None)
def convert_uniprot_to_sequence(uniprot_id) -> str | None:
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


@lru_cache(maxsize=None)
def catalytic_activity(uniprot_id) -> list[str] | None:
    """
    Retrieves the EC (Enzyme Commission) numbers associated with the catalytic activity of a given UniProt ID.

    Parameters:
        uniprot_id (str): The UniProt identifier for the protein of interest.

    Returns:
        list[str] or None: A list of EC numbers if found, otherwise None.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?fields=cc_catalytic_activity"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        ec_numbers = []
        for comment in data.get('comments', []):
            if comment.get('commentType') == 'CATALYTIC ACTIVITY':
                reaction = comment.get('reaction', {})
                ec_number = reaction.get('ecNumber')
                if ec_number:
                    ec_numbers.append(ec_number)
        if len(ec_numbers) != 0:
            return ec_numbers
    else:
        # logging.warning(f"No catalytic activity found for UniProt ID {uniprot_id}")
        return None


def identify_catalytic_enzyme(lst_uniprot_ids, ec) -> str | None:
    """
    Identifies the catalytic enzyme from a list of UniProt IDs for a given EC number.

    Parameters:
        lst_uniprot_ids (str): A semicolon-separated string of UniProt IDs representing enzyme candidates.
        ec (str): The Enzyme Commission (EC) number to match against the catalytic activity.

    Returns:
        str or None: The UniProt ID of the catalytic enzyme if exactly one match is found; 
                     None if no match or multiple matches are found.
    """ 
    enzymes_model = lst_uniprot_ids.split(';')
    catalytic_enzyme = []
    for enzyme in enzymes_model:
        if catalytic_activity(enzyme):
            if ec in catalytic_activity(enzyme):
                catalytic_enzyme.append(enzyme)
    if catalytic_enzyme == []:
        logging.warning(f"{ec}: No catalytic enzyme found for the complex {lst_uniprot_ids}.")
        catalytic_enzyme = None 
    elif len(catalytic_enzyme) > 1:
        logging.warning(f"{ec}: Multiple catalytic enzymes found for the complex {lst_uniprot_ids}.")
        catalytic_enzyme = ';'.join(catalytic_enzyme)
    else:
        catalytic_enzyme = catalytic_enzyme[0]
    return catalytic_enzyme


# if __name__ == "__main__":
    # Test : Send a request to UniProt API
    # uniprot_id = "Q16774"
    # seq = convert_uniprot_to_sequence(uniprot_id)
    # print(seq)

    # Test : Check if catalytic activity is retrieved
    # response = catalytic_activity("P06959")
    # print(response)