# WILDkCAT

**WILDkCAT** is a set of scripts designed to extract, retrieve, and predict enzyme turnover numbers (**kcat**) for genome-scale metabolic models.   

---

## Scripts Overview

### `extract_kcat.py`
- Verifies whether the reaction EC number exists.  
- Retains inputs where reaction-associated genes/enzymes are not supported by KEGG.  
- Retains inputs where no enzymes are provided by the model.  

---

### `retrieve_kcat.py`
- If multiple enzymes are provided, searches UniProt for catalytic activity.  
- If multiple catalytic enzymes are identified, store all.
- When multiple enzymes are found, computes identity percentages relative to the identified catalytic enzyme.  
- Applies Arrhenius correction to values within the appropriate pH range.  
- For rows with multiple scores, selects:
  - The best score  
  - The highest identity percentage  
  - The lowest kcat value  

---

### `predict_kcat.py`
- If multiple enzymes are provided, searches UniProt for catalytic activity.  
- Skips entries missing KEGG compound IDs.  

---

## TODO
- [ ] generate_reports.py: Optimize the functions from generate_reports.py to have html plots
- [ ] catapro.py: Uniprot quieries can be send as batches
- [ ] sabio_rk_api.py: Fix **SABIO-RK** 
- [ ] catapro.py: Move **PubChem API** queries to a dedicated module ? 
- [ ] predict_kcat.py : Integrate **TurNuP** kcat prediction 
- [*] generate_reports.py - general_report(): Add the coverage of the reactions with kcat values in the model

---
