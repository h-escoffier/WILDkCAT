import logging

from wildkcat import run_extraction, run_retrieval, run_prediction_part1, run_prediction_part2, generate_summary_report


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Extraction
    run_extraction("model/e_coli_core.json", "output/ecoli_core/e_coli_core_kcat.tsv")
    # Retrieval
    run_retrieval(
        kcat_file_path="output/ecoli_core/e_coli_core_kcat.tsv",
        output_path="output/ecoli_core/e_coli_core_kcat_retrieved.tsv",
        organism="Escherichia coli",
        temperature_range=(20, 40),
        pH_range=(6.5, 7.5),
        database='brenda'
    ) 
    # Prediction Part 1
    run_prediction_part1("output/ecoli_core/e_coli_core_kcat_retrieved.tsv", 8, "output/ecoli_core/machine_learning/ecoli_catapro_input.csv") 
    
    # run_prediction_part2("output/ecoli_core/e_coli_core_kcat_retrieved.tsv", 
    #                      "output/ecoli_core/machine_learning/ecoli_catapro_output.csv", 
    #                      "output/ecoli_core/machine_learning/ecoli_catapro_input_substrates_to_smiles.tsv", 
    #                      8, 
    #                      "output/ecoli_core/e_coli_core_kcat_full.tsv")
    
    # generate_summary_report("model/e_coli_core.json", "output/ecoli_core/e_coli_core_kcat_full.tsv")
