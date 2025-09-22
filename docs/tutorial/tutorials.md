# Example on _E. coli_ Core Model

This page shows a step-by-step example of the WILDkCAT pipeline on the _E. coli_ core model.

---

## Prerequisites

- Install [WILDkCAT](../installation.md) from PyPI
- Install [CataPro](https://github.com/zchwang/CataPro) to predict kcat values using machine learning
- Download the [E. coli core model](http://bigg.ucsd.edu/static/models/e_coli_core.json)

!!! note 

    All the files used and created in this tutorial are available in the [output folder](https://github.com/h-escoffier/WILDkCAT/tree/main/output) of the WILDkCAT repository

---

## 0 — Create output directories

Ensure all output folders exist to avoid write errors.

```python
import os

os.makedirs("output", exist_ok=True)
os.makedirs("output/machine_learning", exist_ok=True)
```

---

## 1 — Extract kcat values from _E. coli_ core model

First, for each combination of reaction, enzyme, and substrate(s) in the model, create a TSV file. 

Each row corresponds to a unique combination of reaction, enzyme, and substrate(s) and will be used to retrieve experimental kcat values from BRENDA and SABIO-RK in the next step.

```python
from wildkcat import run_extraction

run_extraction(
    model_path="model/e_coli_core.json",
    output_path="output/e_coli_core_kcat.tsv"
)
```

Example of the output file:

| rxn | rxn_kegg | ec_code | direction | substrates_name | substrates_kegg | products_name | products_kegg | genes | uniprot | catalytic_enzyme | warning |
| :-- | :------- | :------ | :-------- | :-------------- | :-------------- | :------------ | :------------ | :---- | :------ | :--------------- | :------ |
| PFK |          | 2.7.1.11 | forward | ATP C10H12N5O13P3;D-Fructose 6-phosphate | C00002;C05345 | ADP C10H12N5O10P2;D-Fructose 1,6-bisphosphate;H+ | C00008;C00354;C00080 | b3916 | P0A796 | P0A796 |  |
| ALCD2x | R00754 | 1.1.1.71 | forward | Ethanol;Nicotinamide adenine dinucleotide | C00469;C00003 | Acetaldehyde;H+;Nicotinamide adenine dinucleotide - reduced | C00084;C00080;C00004 | b0356 | P25437 | P25437 |  |

[View the generated report](extract_ecoli_report.html)

---

## 2 — Retrieve experimental kcat values from BRENDA and SABIO-RK

This function searches for experimentally measured turnover numbers (kcat values) in the BRENDA and/or SABIO-RK databases for the kcats listed in the input file. 
The retrieved values are filtered based on organism, temperature, and pH conditions. The closest matching kcat values are saved to the output file.

```python
from wildkcat import run_retrieval

run_retrieval(
    kcat_file_path="output/e_coli_core_kcat.tsv",
    output_path="output/e_coli_core_kcat_retrieved.tsv",
    organism="Escherichia coli",
    temperature_range=(20, 40),
    pH_range=(6.5, 7.5),
    database='both'
    )
```

Example of the output file:

| rxn | rxn_kegg | ec_code | direction | substrates_name | substrates_kegg | products_name | products_kegg | genes | uniprot | catalytic_enzyme | warning | kcat | matching_score | kcat_substrate | kcat_organism | kcat_enzyme | kcat_temperature | kcat_ph | kcat_variant | kcat_db | kcat_id_percent | kcat_organism_score |
| :-- | :------- | :------ | :-------- | :-------------- | :-------------- | :------------ | :------------ | :---- | :------ | :--------------- | :------ | :--- | :------------- | :------------- | :------------ | :---------- | :--------------- | :------ | :----------- | :------ | :-------------- | :------------------ |
| PFK |          | 2.7.1.11 | forward | ATP C10H12N5O13P3;D-Fructose 6-phosphate | C00002;C05345 | ADP C10H12N5O10P2;D-Fructose 1,6-bisphosphate;H+ | C00008;C00354;C00080 | b3916 | P0A796 | P0A796 | |  0.016 | 1 | D-fructose 6-phosphate | Escherichia coli | P0A796 | 30.0 | 7.2 |  | brenda | 100.0 | 0.0 |
| ALCD2x | R00754 | 1.1.1.71 | forward | Ethanol;Nicotinamide adenine dinucleotide | C00469;C00003 | Acetaldehyde;H+;Nicotinamide adenine dinucleotide - reduced | C00084;C00080;C00004 | b0356 | P25437 |P25437 | | 13.9 | 7 | ethanol | Acinetobacter calcoaceticus |  |  |  |  | brenda |  | 4.0 |

[View the generated report](retrieve_ecoli_report.html)

---

## 3 — Predict missing kcat values using machine learning

### 3.1 - Prepare input file for CataPro

Prepare the input file for CataPro by filtering out the kcat entries that were not found in the previous step and below a limit score (`limit_matching_score`). The resulting file will be used to predict missing kcat values using machine learning.

The function also generates a file named `output_path_substrates_to_smiles.csv` that maps substrate names to their corresponding SMILES. This file will be used to match back the predicted kcat values to the original kcat entries after running CataPro.

```python
from wildkcat import run_prediction_part1

run_prediction_part1(
    kcat_file_path="output/e_coli_core_kcat_retrieved.tsv", 
    output_path="output/machine_learning/ecoli_catapro_input.csv",
    limit_matching_score=6
    )
```

The output file `ecoli_catapro_input.csv` is formatted according to the requirements of [CataPro](https://github.com/zchwang/CataPro), meaning it can be directly used as input for kcat prediction.

!!! note 

    Before running predictions, make sure you have installed CataPro by following the installation instructions provided in their [GitHub repository](https://github.com/zchwang/CataPro?tab=readme-ov-file#create-the-catapro-environment).

Once installed, you can run CataPro with the following command:

```bash 
python predict.py \
        -inp_fpath output/machine_learning/ecoli_catapro_input.csv \
        -model_dpath models \
        -batch_size 64 \
        -device cuda:0 \
        -out_fpath ecoli_catapro_output.csv
```

[View the generated report](predict_ecoli_report.html)


### 3.2 - Integrate CataPro predictions

After running CataPro with the prepared input file, integrate the predicted kcat values back into the original kcat entries. The function matches the predicted values to the original entries using the substrate names and SMILES mapping file generated in the previous step.

```python
from wildkcat import run_prediction_part2

run_prediction_part2(
    kcat_file_path="output/e_coli_core_kcat_retrieved.tsv", 
    catapro_predictions_path="output/machine_learning/ecoli_catapro_output.csv", 
    substrates_to_smiles_path="output/machine_learning/ecoli_catapro_input_substrates_to_smiles.tsv", 
    output_path="output/e_coli_core_kcat_full.tsv",
    limit_matching_score=6
    )
```

Example of the output file:

| rxn | rxn_kegg | ec_code  | direction | substrates_name | substrates_kegg  | products_name | products_kegg | genes | uniprot | catalytic_enzyme | warning | kcat | db | matching_score | kcat_substrate | kcat_organism | kcat_enzyme | kcat_temperature | kcat_ph | kcat_variant | kcat_id_percent |
| :-- | :------- | :------- | :-------- | :-------------- | :--------------- | :------------ | :-------------| :---- | :------ | :--------------- | :------ | :--- | :- | :------------- | :------------- | :------------ | :---------- | :--------------- | :------ | :----------- | :-------------- |
| PFK |          | 2.7.1.11 | forward   | ATP C10H12N5O13P3; D-Fructose 6-phosphate | C00002; C05345 | ADP C10H12N5O10P2; D-Fructose 1,6-bisphosphate; H+ | C00008; C00354; C00080 | b3916 | P0A796 | P0A796 | | 0.016 | brenda  | 1 | D-fructose 6-phosphate | Escherichia coli | P0A796 | 30.0 | 7.2 |  | 100.0 |
| ALCD2x | R00754 | 1.1.1.71 | forward | Ethanol;Nicotinamide adenine dinucleotide | C00469;C00003 | Acetaldehyde;H+;Nicotinamide adenine dinucleotide - reduced | C00084;C00080;C00004 | b0356 | P25437 | P25437| | 16.0905 | catapro |  |  |  |  |  |  |  |  |

---

## 4 — Generate summary report

The final output file `e_coli_core_kcat_full.tsv` contains both experimentally retrieved and machine learning predicted kcat values for each combination of reaction, enzyme, and substrate(s) in the _E. coli_ core model. This file can be used for integration into enzyme-constrained metabolic models.

The result can be visualized and summarized using the function `generate_summary_report`: 

```python
from wildkcat.visualization import generate_summary_report

generate_summary_report(
    model_path="model/e_coli_core.json", 
    kcat_file_path="output/e_coli_core_kcat_full.tsv"
    )
```

[View the generated report](general_ecoli_report.html)
