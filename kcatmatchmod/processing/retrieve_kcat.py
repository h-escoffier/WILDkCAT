import logging
from api import run_sabio_rk

def retrieve_kcat_from_brenda(): 
    pass 


def retrieve_kcat_from_sabio_rk(kcat_path, organism, temperature_range, pH_range):
    run_sabio_rk(kcat_path,
                 organism,
                 temperature_range,
                 pH_range)


def retrieve_kcat_from_uniprot(): 
    pass 


def retrieve_kcat_from_catapro(): 
    pass 


def retrieve_kcat_from_turnip():
    pass


def run_retrieve_kcat(kcat_path, organism, temperature_range, pH_range):
    # Retrieve kcat values from SABIO-RK
    df = retrieve_kcat_from_sabio_rk(kcat_path, organism, temperature_range, pH_range)
    # Retrieve kcat values from BRENDA
    df = retrieve_kcat_from_brenda()
    # Retrieve kcat values from UniProt
    df = retrieve_kcat_from_uniprot()
    # Retrieve kcat values from CaTaPro
    df = retrieve_kcat_from_catapro()
    # Save the output DataFrame to a file
    df.to_csv(kcat_path, sep="\t", index=False)
