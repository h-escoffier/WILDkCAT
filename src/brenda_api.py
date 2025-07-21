from zeep import Client, Settings
from dotenv import load_dotenv
import hashlib
import os 


load_dotenv() 


def get_turnover_number_brenda(
    ec_number: str = "",
    organism: str = "",
    turnover_number: str = "",
    turnover_number_maximum: str = "",
    substrate: str = "",
    commentary: str = "",
    ligand_structure_id: str = "",
    literature: str = ""
):
    """
    Queries the BRENDA SOAP API to retrieve turnover number values for a given enzyme.

    Args:
        email (str): Your registered email address with BRENDA.
        password (str): Your BRENDA password or API key.
        ec_number (str, optional): Enzyme Commission (EC) number (e.g., '1.1.1.1').
        organism (str, optional): Organism name (e.g., 'Homo sapiens').
        turnover_number (str, optional): Specific turnover number to filter by.
        turnover_number_maximum (str, optional): Maximum turnover number to filter by.
        substrate (str, optional): Substrate involved.
        commentary (str, optional): Commentary associated with the entry.
        ligand_structure_id (str, optional): Ligand structure ID.
        literature (str, optional): Literature reference(s), comma-separated if multiple.

    Returns:
        list[dict]: A list of dictionaries containing turnover number entries.

    Raises:
        ValueError: If neither `ec_number` nor `organism` is specified.
    """

    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")

    if not ec_number and not organism:
        raise ValueError("You must specify at least one of 'ec_number' or 'organism'.")

    # Call the SOAP API
    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    settings = Settings(strict=False)

    parameters = [
        email,
        password,
        f"ecNumber*{ec_number}",
        f"organism*{organism}",
        f"turnoverNumber*{turnover_number}",
        f"turnoverNumberMaximum*{turnover_number_maximum}",
        f"substrate*{substrate}",
        f"commentary*{commentary}",
        f"ligandStructureId*{ligand_structure_id}",
        f"literature*{literature}",
    ]

    client = Client(wsdl, settings=settings)
    result = client.service.getTurnoverNumber(*parameters)

    return result


# Notes : 
# The BRENDA API gives multiple entries for the same EC number and substrate.
# Remove the co-substrates

