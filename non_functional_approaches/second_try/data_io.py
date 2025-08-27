import pandas as pd 


def format_file_for_catapro(kcat_path, output_file):
    """
    Format an enzyme file by extracting and renaming specific columns.
    
    Parameters:
    - input_file (str): Path to the original TSV file.
    - output_file (str): Path to save the formatted TSV file.
    """
    # Load the original file
    df = pd.read_csv(kcat_path, sep='\t')
    # test 
    df = df.head(10)
    # Create a new DataFrame with required columns
    formatted_df = pd.DataFrame({
        'Enzyme_id': df['Uniprot'],
        'type': 'wild',
        'sequence': df['Sequence'],
        'smiles': df['SMILES']
    })
    formatted_df.to_csv(output_file, sep=',', index=True)


if __name__ == "__main__":
    format_file_for_catapro("output/CataPro/kcat_Human-GEM_smiles_sequences.tsv", "output.tsv")