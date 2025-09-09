import logging

from wildkcat import run_extraction, run_retrieval, run_prediction_part1, run_prediction_part2, generate_summary_report


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Extraction
    run_extraction(
        model_path="model/e_coli_core.json", 
        output_path="output/ecoli_core/e_coli_core_kcat.tsv"
    )
    
    # Retrieval
    run_retrieval(
        kcat_file_path="output/ecoli_core/e_coli_core_kcat.tsv",
        output_path="output/ecoli_core/e_coli_core_kcat_retrieved.tsv",
        organism="Escherichia coli",
        temperature_range=(20, 40),
        pH_range=(6.5, 7.5),
        database='both'
    ) 

    # # Prediction Part 1
    run_prediction_part1(
        kcat_file_path="output/ecoli_core/e_coli_core_kcat_retrieved.tsv", 
        output_path="output/ecoli_core/machine_learning/ecoli_catapro_input.csv",
        limit_matching_score=8 
    )

    # # Prediction Part 2
    run_prediction_part2(
        kcat_file_path="output/ecoli_core/e_coli_core_kcat_retrieved.tsv", 
        catapro_predictions_path="output/ecoli_core/machine_learning/ecoli_catapro_output.csv", 
        substrates_to_smiles_path="output/ecoli_core/machine_learning/ecoli_catapro_input_substrates_to_smiles.tsv", 
        output_path="output/ecoli_core/e_coli_core_kcat_full.tsv",
        limit_matching_score=8
    )

    # # Summary Report
    generate_summary_report(
        model_path="model/e_coli_core.json", 
        kcat_file_path="output/ecoli_core/e_coli_core_kcat_full.tsv"
    )
