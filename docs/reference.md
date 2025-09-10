# Main functions 

The **WILDkCAT** package is organized into modules: 

1. **Extraction** : Extraction of *kcat* values from the provided model

2. **Retrieval** : Retrieval of *kcat* values using curated databases (BRENDA and SABIO-RK)

3. **Prediction** : Prediction of missing and low confidence *kcat* values using ML-based CataPro model

4. **Summary** : Generates an HTML report summarizing the percentage and quality of kcat values identified for the model, along with their data sources.

---

## Extraction 
::: wildkcat.processing.extract_kcat.run_extraction

---

## Retrieval
:::wildkcat.processing.retrieve_kcat.run_retrieval

---

## Prediction
:::wildkcat.processing.predict_kcat.run_prediction_part1
:::wildkcat.processing.predict_kcat.run_prediction_part2

---

## Summary report
:::wildkcat.processing.summary.generate_summary_report

---

## Matching process and scoring

The matching process is designed to select the most appropriate *kcat* value when multiple candidates are available.  
Each candidate is first assigned a score based on several criteria, such as:  

- **kcat specific criteria:**
    - Substrate
    - Catalytic enzyme(s)
- **General criteria:**
    - Organism
    - Temperature
    - pH

If two or more candidates receive the same score, tie-breaking rules are applied in the following order:  

1. **Enzyme sequence identity** – the value associated with the most similar protein sequence is preferred.  
2. **Organism proximity** – preference is given to *kcat* values measured in organisms closest to the target species.  
3. **Minimal *kcat* value** – if ambiguity remains, the smallest *kcat* value is chosen. 

:::wildkcat.utils.matching
heading_level: 3
