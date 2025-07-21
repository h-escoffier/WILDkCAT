import requests
import pandas as pd 
from io import StringIO


def get_turnover_number_sabio(
        ec_number: str = "",
        organism: str = "", 
        kegg_reaction_id: str = ""
):
    """
    TODO 
    """
    base_url = 'https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'
    parameters = 'https://sabiork.h-its.org/entry/exportToExcelCustomizable'
    entryIDs = []


    # ask SABIO-RK for all EntryIDs matching a query
    query_parts = ['Parametertype:"kcat"']
    # query_parts = []
    if ec_number:
        query_parts.append(f'ECNumber:"{ec_number}"')
    if organism:
        query_parts.append(f'Organism:"{organism}"')
    if kegg_reaction_id:
        query_parts.append(f'KeggReactionID:"{kegg_reaction_id}"')
    query_string = ' AND '.join(query_parts)
    query = {'format':'txt', 'q':query_string}

    # make GET request
    request = requests.get(base_url, params = query)
    request.raise_for_status() # raise if 404 error

    # each entry is reported on a new line
    entryIDs = [int(x) for x in request.text.strip().split('\n')]
    print('%d matching entries found.' % len(entryIDs))

    # encode next request, for parameter data given entry IDs
    data_field = {'entryIDs[]': entryIDs}
    query = {'format':'tsv', 'fields[]':['EntryID', 'Organism', 'UniprotID','ECNumber', 'Parameter', 'ReactomeReactionID']}
    
    print(request.url)

    # make POST request
    request = requests.post(parameters, params=query, data=data_field)
    request.raise_for_status()

    exit()

    # Load response into a pandas DataFrame
    df = pd.read_csv(StringIO(request.text), sep='\t')

    # Filter for only 'kcat' parameter types
    kcat_df = df[df['parameter.name'].str.lower() == 'kcat'].reset_index(drop=True)
    return kcat_df


# if __name__ == "__main__":
#     # Example usage
#     ec_number = "1.1.1.1"
#     organism = "Homo sapiens"
#     kegg_reaction_id = "R00754"
#     output = get_turnover_number_sabio(ec_number, organism, kegg_reaction_id)
#     output.to_csv('test.tsv', sep="\t", index=False)

import requests

# First: get entry IDs
response = requests.get(
    'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs',
    params={"q": 'ECNumber:"2.7.1.11" AND Organism:"Escherichia coli"',
        "format": 'txt'}
)

response.raise_for_status()
entryIDs = [int(x) for x in response.text.strip().split('\n')]
print('%d matching entries found.' % len(entryIDs))

# Second: download the data for those IDs
# encode next request, for parameter data given entry IDs
data_field = {'entryIDs[]': entryIDs}
query = {'format':'tsv', 
        #  'fields[]':['EntryID', 'Organism', 'UniprotID','ECNumber', 'Parameter', 'ReactomeReactionID']}
        'fields[]':['Parametertype', 'DateSubmitted','PubMedID', 'Parameter']}
# make POST request
request = requests.post('https://sabiork.h-its.org/entry/exportToExcelCustomizable', params=query, data=data_field)
request.raise_for_status()

print(request.text)
