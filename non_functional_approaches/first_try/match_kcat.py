import pandas as pd
import ast

# Note : Matching kcat values for Human-GEM gives no results. 

def match_kcat(hgem_file, kcat_file, output_file=None):
    # Load the data
    hgem_df = pd.read_csv(hgem_file, sep='\t')
    kcat_df = pd.read_csv(kcat_file, sep='\t')

    # Function to safely parse CHEBI lists
    def parse_chebi_list(x):
        if isinstance(x, str) and x.startswith('['):
            return set(ast.literal_eval(x))
        elif isinstance(x, str) and ';' in x:
            return set(x.split(';'))
        elif isinstance(x, str):
            return {x}
        return set()

    # Preprocess the kcat data
    kcat_df['Uniprot ID'] = kcat_df['Uniprot ID'].astype(str)
    kcat_df['substrates'] = kcat_df['substrate CHEBI IDs'].apply(parse_chebi_list)
    kcat_df['products'] = kcat_df['product CHEBI IDs'].apply(parse_chebi_list)

    # Preprocess the Human-GEM data
    hgem_df['Uniprot'] = hgem_df['Uniprot'].astype(str)
    hgem_df['substrates'] = hgem_df['chebi_Substrate'].apply(parse_chebi_list)
    hgem_df['products'] = hgem_df['chebi_Product'].apply(parse_chebi_list)

    # Function to find a matching kcat
    def find_kcat(row):
        matches = kcat_df[
            (kcat_df['Uniprot ID'] == row['Uniprot']) &
            (kcat_df['substrates'] == row['substrates']) &
            (kcat_df['products'] == row['products'])
        ]
        if not matches.empty:
            return matches.iloc[0]['kcat [1/sec]']
        return None

    # Apply matching function
    hgem_df['kcat [1/sec]'] = hgem_df.apply(find_kcat, axis=1)

    # Optionally write the output
    if output_file:
        hgem_df.to_csv(output_file, sep='\t', index=False)

    return hgem_df


# if __name__ == "__main__":
#     print("start")
#     result_df = match_kcat("output/Human-GEM_kcat_Chebi_cleaned.tsv", "data/kcat/Uniprot_kcat.tsv", output_file="output/Human-GEM_Chebi_with_kcat.tsv")
#     print("end")