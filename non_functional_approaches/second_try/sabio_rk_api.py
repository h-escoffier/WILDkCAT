import requests
import pandas as pd 
from io import StringIO
from tqdm import tqdm

def get_turnover_number_sabio(
        ec_number: str = "",
        organism: str = "", 
        kegg_reaction_id: str = "",
        enzyme_variant: str = "", 
        ph_range: tuple = (7.0, 7.4),
        temperature_range: tuple = (37.0, 37.0)
):
    """
    TODO 
    """
    base_url = 'https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'
    parameters = 'https://sabiork.h-its.org/entry/exportToExcelCustomizable'
    entryIDs = []

    # 1. Design query and retrieve entry IDs
    query_parts = ['Parametertype:"kcat"']
    if ec_number:
        query_parts.append(f'ECNumber:"{ec_number}"')
    if organism:
        query_parts.append(f'Organism:"{organism}"')
    if kegg_reaction_id:
        query_parts.append(f'KeggReactionID:"{kegg_reaction_id}"')
    query_string = ' AND '.join(query_parts)
    query = {'format':'txt', 'q':query_string}

    # Make GET request
    request = requests.get(base_url, params=query)
    request.raise_for_status()  # Raise if 404 error
    if request.text == "no data found":
        # print('%s: No data found for the query.' % query_string)  # TODO : add logging or printlevel arg
        return pd.DataFrame()  # Return empty DataFrame if no data found

    entryIDs = [int(x) for x in request.text.strip().split('\n')]
    # print('%d matching entries found.' % len(entryIDs)) # TODO : add logging or printlevel arg

    # 2. Design the request to retrieve informations matching the entryIDs
    data_field = {'entryIDs[]': entryIDs}
    # Possible fields to retrieve:
    # EntryID, Reaction, Buffer, ECNumber, CellularLocation, UniProtKB_AC, Tissue, Enzyme Variant, Enzymename, Organism
    # Temperature, pH, Activator, Cofactor, Inhibitor, KeggReactionID, KineticMechanismType, Other Modifier, Parameter,
    # Pathway, Product, PubMedID, Publication, Rate Equation, SabioReactionID, Substrate
    query = {'format':'tsv', 'fields[]':['EntryID', 'Reaction', 'Buffer', 'ECNumber', 'CellularLocation', 'UniProtKB_AC', 'Tissue', 'Enzyme Variant', 'Enzymename', 'Organism',
                                         'Temperature', 'pH', 'Activator', 'Cofactor', 'Inhibitor', 'KeggReactionID', 'KineticMechanismType', 'Parameter', 'Product', 'Substrate']}

    # Make POST request
    request = requests.post(parameters, params=query, data=data_field)
    request.raise_for_status()

    # format the response 
    df = pd.read_csv(StringIO(request.text), sep='\t')

    # print(df['Temperature'])
    # print(df['pH'])

    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['pH'] = pd.to_numeric(df['pH'], errors='coerce')

    # print(df['Temperature'])
    # print(df['pH'])

    # Filters
    df = df[df['parameter.name'].str.lower() == 'kcat'].reset_index(drop=True) # Keep only kcat parameters
    if enzyme_variant:
        df = df[df['Enzyme Variant'].str.lower() == enzyme_variant].reset_index(drop=True) # Keep only the specified enzyme variant
    df = df[(df['Temperature'] >= temperature_range[0]) & (df['Temperature'] <= temperature_range[1])].reset_index(drop=True) # Keep only the specified temperature
    df = df[(df['pH'] >= ph_range[0]) & (df['pH'] <= ph_range[1])].reset_index(drop=True) # Keep only the specified pH range
    return df


def use_sabio_rk_api(rxns_path, organism, ph_range, temperature_range, enzyme_variant, output_path=None):
    """
    Reads a TSV file with columns 'ec_code' and/or 'kegg_reaction_id' and queries SABIO-RK for kcat values.
    Returns a DataFrame where each reaction can have multiple rows (one per SABIO entry).
    """
    rxns_df = pd.read_csv(rxns_path, sep='\t', dtype=str)

    results = []

    for _, row in tqdm(iterable=rxns_df.iterrows(), total=rxns_df.shape[0], desc="Querying SABIO-RK"):
        ec = row.get('ec_code', '') or ''
        kegg = row.get('kegg_reaction_id', '') or ''

        print(f"Processing EC: {ec}, KEGG: {kegg}")

        # Call SABIO API for this reaction
        df_kcat = get_turnover_number_sabio(
            ec_number=ec,
            organism=organism,
            kegg_reaction_id=kegg,
            enzyme_variant=enzyme_variant,
            ph_range=ph_range,
            temperature_range=temperature_range
        )

        print(f"Found {df_kcat.shape[0]} entries for EC: {ec}, KEGG: {kegg}")

        if not df_kcat.empty:
            # Add identifiers from the input file to keep mapping
            df_kcat.insert(0, 'ec_code', ec)
            df_kcat.insert(1, 'kegg_reaction_id', kegg)
            results.append(df_kcat)

    # Concatenate all results (each SABIO entry = one row)
    if results:
        final_df = pd.concat(results, ignore_index=True)
    else:
        # Return empty DataFrame with consistent columns
        final_df = pd.DataFrame(columns=['ec_code', 'kegg_reaction_id', 'EntryID', 'Reaction', 'Buffer', 'ECNumber', 'CellularLocation', 'UniProtKB_AC', 'Tissue', 'Enzyme Variant', 'Enzymename', 'Organism',
                                         'Temperature', 'pH', 'Activator', 'Cofactor', 'Inhibitor', 'KeggReactionID', 'KineticMechanismType', 'Other Modifier', 'Parameter',
                                         'Pathway', 'Product', 'PubMedID', 'Publication', 'Rate Equation', 'SabioReactionID', 'Substrate'])
    if output_path:
        final_df.to_csv(output_path, sep='\t', index=False)
        # print(f"Results saved to {output_path}") TODO : add logging or printlevel arg
    
    return final_df

    
if __name__ == "__main__":
    # Example usage
    df = use_sabio_rk_api(
        rxns_path='output/EColiCore/sabio_rk.tsv',
        organism='Escherichia coli',
        ph_range=(7.0, 8.5),
        temperature_range=(21.0, 37.0),
        enzyme_variant='wildtype',
        output_path='output/EColiCore/sabio_rk_results.tsv'
    )