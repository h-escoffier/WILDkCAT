from wildkcat import run_extraction, run_retrieval, run_prediction_part1, run_prediction_part2, generate_summary_report


if __name__ == "__main__":
    # Extraction
    run_extraction(
        model_path="model/e_coli_core.json", 
        output_folder="output"
    )
    
    # Retrieval
    run_retrieval(
        output_folder="output",
        organism="Escherichia coli",
        temperature_range=(20, 45),
        pH_range=(4, 8),
        database="both"
    ) 

    # Prediction (OPTIONAL)
    ## Prediction Part 1
    run_prediction_part1(
        output_folder="output",
        limit_penalty_score=9
    )

    ## Prediction Part 2
    run_prediction_part2(
        output_folder="output", 
        catapro_predictions_path="output/machine_learning/catapro_output.csv", 
        limit_penalty_score=9
    )

    # Summary Report
    generate_summary_report(
        model_path="model/e_coli_core.json", 
        output_folder="output"
    )
