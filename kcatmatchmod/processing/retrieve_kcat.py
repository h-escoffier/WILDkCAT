import logging

from kcatmatchmod.api.sabio_rk_api import run_sabio_rk
from kcatmatchmod.api.brenda_api import run_brenda


def retrieve_kcat_from_brenda(kcat_path, organism, temperature_range, pH_range):
    # kcat_file_path, output_path, organism, temperature_range, pH_range, variant = "wildtype", report=True
    df = run_brenda(kcat_path,
                    kcat_path.replace('.tsv', '_brenda.tsv'),
                    organism,
                    temperature_range,
                    pH_range
                    )
    return df


def retrieve_kcat_from_sabio_rk(kcat_path, organism, temperature_range, pH_range):
    # kcat_file_path, output_path, organism, temperature_range, pH_range, variant = "wildtype", report=True
    df = run_sabio_rk(kcat_path,
                      kcat_path.replace('.tsv', '_sabio.tsv'),
                      organism,
                      temperature_range,
                      pH_range
                     )
    return df


def retrieve_kcat_from_uniprot(): 
    pass 


def retrieve_kcat_from_catapro(): 
    pass 


def retrieve_kcat_from_turnip():
    pass


def run_retrieve_kcat(kcat_path, organism, temperature_range, pH_range):
    # Retrieve kcat values from SABIO-RK
    logging.info("Retrieving kcat values from SABIO-RK")
    df = retrieve_kcat_from_sabio_rk(kcat_path, organism, temperature_range, pH_range)
    # Retrieve kcat values from BRENDA
    logging.info("Retrieving kcat values from BRENDA")
    df = retrieve_kcat_from_brenda(kcat_path, organism, temperature_range, pH_range)
    # Retrieve kcat values from UniProt
    df = retrieve_kcat_from_uniprot()
    # Retrieve kcat values from CaTaPro
    df = retrieve_kcat_from_catapro()
    # Save the output DataFrame to a file
    df.to_csv(kcat_path, sep="\t", index=False)
