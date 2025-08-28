import typer
from dotenv import load_dotenv
from wildkcat import run_extraction, run_retrieval, run_prediction_part1, run_prediction_part2, generate_summary_report


load_dotenv()


app = typer.Typer(help="WILDkCAT CLI - Extract, Retrieve and Predict kcat values for a metabolic model.")


@app.command()
def extraction(
    model_path: str, 
    output_path: str, 
    report: bool = True
):
    """
    Extracts kcat-related data from a metabolic model and generates output files and an optional HTML report.

    Parameters:

        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).

        output_path (str): Path to the output file (TSV format).

        report (bool, optional): Whether to generate an HTML report (default: True).
    """
    run_extraction(model_path=model_path, output_path=output_path, report=report)
    typer.echo(f"Extraction finished. Output saved at {output_path}")


@app.command()
def retrieval(
    kcat_file_path: str,
    output_path: str,
    organism: str,
    temperature_range: tuple[float, float],
    ph_range: tuple[float, float],
    database: str = 'both',
    report: bool = True
):
    """
    Retrieves closests kcat values from specified databases for entries in a kcat file, applies filtering criteria, 
    and saves the results to an output file.
    
    Parameters:
    
        kcat_file_path (str): Path to the input kcat file.

        output_path (str): Path to save the output file with retrieved kcat values.

        organism (str): Organism name.

        temperature_range (tuple): Acceptable temperature range for filtering (min, max).

        ph_range (tuple): Acceptable pH range for filtering (min, max).

        database (str, optional): Specifies which database(s) to query for kcat values. Options are 'both' (default), 'brenda', or 'sabio_rk'.

        report (bool, optional): Whether to generate an HTML report using the retrieved data (default: True).       
    """
    run_retrieval(
        kcat_file_path=kcat_file_path,
        output_path=output_path,
        organism=organism,
        temperature_range=temperature_range,
        pH_range=ph_range,
        database=database,
        report=report
    )
    typer.echo(f"Retrieval finished. Output saved at {output_path}")


@app.command()
def prediction_part1(
    kcat_file_path: str, 
    output_path: str, 
    limit_matching_score: int,
    report: bool = True
):
    """
    Processes kcat data file to generate input files for CataPro prediction.
    Optionally, it can produce a summary report of the processed data.

    Parameters:
    
        kcat_file_path (str): Path to the input kcat data file.

        output_path (str): Path to save the generated CataPro input CSV file.

        limit_matching_score (int): Threshold for filtering entries based on matching score.

        report (bool, optional): Whether to generate a report using the retrieved data (default: True). 
    """
    run_prediction_part1(
        kcat_file_path=kcat_file_path,
        output_path=output_path,
        limit_matching_score=limit_matching_score, 
        report=report
    )
    typer.echo(f"Prediction Part 1 finished. Output saved at {output_path}")


@app.command()
def prediction_part2(
    kcat_file_path: str,
    catapro_predictions_path: str,
    substrates_to_smiles_path: str,
    output_path: str,
    limit_matching_score: int
):
    """
    Runs the second part of the kcat prediction pipeline by integrating Catapro predictions,
    mapping substrates to SMILES, formatting the output, and optionally generating a report.
    
    Parameters:
        kcat_file_path (str): Path to the input kcat TSV file.

        catapro_predictions_path (str): Path to the Catapro predictions CSV file.

        substrates_to_smiles_path (str): Path to the TSV file mapping substrates to SMILES.

        output_path (str): Path to save the formatted output TSV file.

        limit_matching_score (float): Threshold for taking predictions over retrieved values.

        report (bool, optional): If True, generates a report (default: True). 
    """
    run_prediction_part2(
        kcat_file_path=kcat_file_path,
        catapro_predictions_path=catapro_predictions_path,
        substrates_to_smiles_path=substrates_to_smiles_path,
        output_path=output_path,
        limit_matching_score=limit_matching_score
    )
    typer.echo(f"Prediction Part 2 finished. Output saved at {output_path}")


@app.command()
def report(model_path: str, kcat_file_path: str):
    """
    Generate a HTML report summarizing the kcat extraction, retrieval and prediction for a given model. 

    Parameters:

        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).
        
        kcat_file_path (str): Path to the final kcat TSV file.
    """
    generate_summary_report(
        model_path=model_path,
        kcat_file_path=kcat_file_path
    )
    typer.echo(f"Summary report generated. Output saved at reports/general_report.html")


if __name__ == "__main__":
    app()
