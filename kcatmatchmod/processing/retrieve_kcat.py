from api import run_sabio_rk

def retrieve_kcat_from_brenda(): 
    pass 


def retrieve_kcat_from_sabio_rk():
    run_sabio_rk(input_file='output/ecoli_kcat.tsv',
                 organism='Escherichia coli',
                 temperature_range=(20, 37),
                 pH_range=(6, 8))


def retrieve_kcat_from_catapro(): 
    pass 


def retrieve_kcat_from_turnip():
    pass


def run_retrieve_kcat(kcat_path, organism, temperature_range, pH_range):
    pass 
