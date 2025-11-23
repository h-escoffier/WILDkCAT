# Explanation of the WILDkCAT Pipeline

## Workflow Overview

![WILDkCAT Workflow](workflow.svg)

---

## 1 - Extract kcat values from a metabolic model

The first step involves extracting potential kcat values from a given metabolic model. 

The output is a TSV file where each row corresponds to a unique combination of reaction, enzyme, and substrate(s). This file serves as the basis for retrieving experimental kcat values from BRENDA and SABIO-RK databases.

### Details of the extraction process

During this step, several verification procedures are carried out to ensure consistency between the model annotations and external databases:

1. **Fields checked during extraction**

    * Reaction identifiers: `model.reaction[i].annotation.get('kegg.reaction')`
    * EC numbers: `model.reaction[i].annotation.get('ec-code')`
    * Compounds: `model.metabolites[i].annotation.get('kegg.compound')`
    * Genes: `model.genes[i].annotation.get('uniprot')`

2. **Verification of EC numbers via KEGG**

    Each EC number associated with a reaction is checked using the KEGG API. This ensures that the EC number is still valid, or if it has been transferred. 

3. **Handling enzyme complexes (multiple enzymes per reaction)**

    * When a reaction is associated with multiple enzymes (enzyme complex), the UniProt API is queried to identify which subunit is the catalytic enzyme.

    * Two dedicated columns (`warning_ec` and `warning_enz`) are used to flag potential issues:
        - `warning_ec`:
            - `missing`: no EC number found for the reaction.
            - `invalid`: EC number is invalid.
            - `transferred`: EC number has been transferred to a new one.
            - (empty): EC number is valid.
        - `warning_enz`:
            - `no_gpr`: no gene-protein-reaction (GPR) association found for the reaction.
            - `none`: no catalytic enzyme identified.
            - `multiple`: more than one possible catalytic enzyme identified.
            - (empty): one clear catalytic enzyme has been identified.

4. **Retrieval of kcat with incomplete information**    

    When an EC number is incomplete or has been transferred (c.f. `warning_ec` column), WILDkCAT performs the retrieval based only on the available enzyme information. These situations reduce the likelihood of finding alternative kcat values.

5. **Removing rows with insufficient information**

    If both the EC number and the catalytic enzyme information are missing for a given reaction, the corresponding row is removed due to insufficient information to assign a kcat value, theses rows are notified in the table generated in the `extract_report.html`.

---

## 2 - Retrieve experimental kcat values from BRENDA and/or SABIO-RK

In this step, experimental kcat values are retrieved from the BRENDA and SABIO-RK databases.

If no exact value is available, the pipeline assigns the closest possible match using a penalty-based scoring system.

### Matching score

The matching score evaluates how well a candidate kcat entry fits the query enzyme and conditions.

* A lower score indicates a better match.
* `0` = best possible match (perfect fit).
* `17` = no reliable match.

### Penalty-based scoring system 

!!! warning 

    Penalty values have been set in collaboration with an expert ([Carole Linster](https://www.uni.lu/lcsb-en/people/carole-linster/)), they can be adjusted based on user feedback. To modify these values, please refer to the `wildkcat/utils/matching.py` file or contact us by creating an issue on [GitHub](https://github.com/h-escoffier/WILDkCAT/issues/new). 

Each criterion adds a penalty if the candidate entry deviates from the query:

| Criterion            | Description                                           | Penalty |
|----------------------|-------------------------------------------------------|---------|
| **Organism**         | Same organism                                         | 0       |
|                      | Different or unknown                                  | 2       |
| **Substrate**        | Same substrate                                        | 0       |
|                      | Different or unknown                                  | 4       |
| **Catalytic enzyme** | Matches expected enzyme                               | 0       |
|                      | Does not match or unknown                             | 3       |
| **Variant**          | Wild-type                                             | 0       |
|                      | Unknown                                               | 1       |
|                      | Mutant                                                | 16      |
| **Temperature**      | Within specified range                                | 0       |
|                      | Corrected via Arrhenius equation                      | 0       |
|                      | Unknown                                               | 1       |
|                      | Outside specified range                               | 2       |
| **pH**               | Within specified range                                | 0       |
|                      | Unknown                                               | 1       |
|                      | Outside specified range                               | 2       |
    
### Details of the retrieval process

1. **Database querying**

    The BRENDA and/or SABIO-RK databases are queried using their respective APIs. Queries are performed based on EC numbers.

2. **Tie-Breaking Strategy**

    When multiple kcat entries share the same matching score, the following criteria are applied sequentially to select the most appropriate value:

    * *Sequence Identity*: Align the enzyme sequences using [BioPython](https://biopython.org)’s `Align.PairwiseAligner()` and prioritize higher sequence identity.
    
    * *Organism Proximity*: Use `Entrez.efetch()` to determine the organism closest to the target species.

    * *Maximal kcat Selection*: If ties remain, select the entry with the highest kcat value, thus avoiding excessively stringent constraints on the model.

3. **Correction of the kcat value using the Arrhenius equation**

    If the temperature at which kcat was measured is outside the desired range and at least two kcat measurements are available (to estimate \( E_a \)), the kcat value can be adjusted using the [Arrhenius equation](https://en.wikipedia.org/wiki/Arrhenius_equation):

    $$
    kcat_{\text{opt}} = kcat_{\text{db}} \times \exp\left(-\frac{E_a}{R} \left(\frac{1}{T_{\text{db}}} - \frac{1}{T_{\text{opt}}}\right)\right)
    $$

    Where:

    - \( E_a \) = Activation energy
    - \( R \) = Universal gas constant (8.314 J/(mol·K))
    - \( T_{db} \) = Temperature at which kcat was measured (in Kelvin)
    - \( T_{opt} \) = Temperature range midpoint (in Kelvin)

    Notes: 

    *  The target temperature \( T_{opt} \) corresponds to the midpoint of the desired optimal range.

4. **Merging rows with same reaction-enzyme combination**

    During the extraction step, multiple rows may be generated if a combination of reaction and enzyme is associated with multiple EC numbers. In such cases, after the retrieval step, these rows are merged to retain only the best kcat value. The column `ec_codes` lists all the EC numbers associated with the reaction-enzyme pair. The column `ec_code` indicates the EC number corresponding to the selected kcat value.


5. **Caching intermediate retrieval results** 

    During the retrieval step, partial results are periodically stored (every 200 rows) in a temporary folder named `cache_retrieval`. This caching mechanism allows the process to resume seamlessly after an interruption (e.g. API failure or network issue) without restarting from the beginning. Each row in the dataset is marked with a processed flag to indicate whether it has already been queried, ensuring that previously completed entries are skipped upon resumption. Once the retrieval completes successfully and the final file is generated, the cache folder is keeped and the `processed` column is removed to keep the output file clean.

---

## 3 - Predict missing kcat values using machine learning

!!! note 

    The prediction step is optional and can be skipped if only experimental kcat values are desired.

When no suitable experimental kcat value is found, the pipeline allows the prediction of kcat values using the CataPro machine learning model.
The predictions rely on reaction substrates (SMILES) and enzyme sequences. 

Predicted kcat values are used to replace experimental values that fall below a threshold, defined by the `limit_matching_score` argument.
If a kcat value cannot be predicted, the best available experimental kcat value is retained.

If multiple substrates are involved in the reaction, the prediction is performed for each substrate, and the lowest predicted kcat value is retained.

### Details of the prediction process

1. **Identification of the enzyme sequence**

    The UniProt API is queried to retrieve the amino acid sequence of the catalytic enzyme from the UniProt ID.

2. **Cofactor Identification**

    The BRENDA API is queried to identify cofactors associated with the reaction. If no cofactors are found, the kcat prediction is performed for all the metabolites involved in the reaction, then the lowest predicted kcat value is retained.

3. **Substrate SMILES Retrieval**

    The PubChem API is queried to retrieve the SMILES representation of each substrate from their KEGG IDs.

---