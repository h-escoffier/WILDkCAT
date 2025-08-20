<p align="center">
  <img src="docs/WILDKCAT.tif" alt="WildKcat Logo"/>
</p>

---

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
- [ ] Fix **SABIO-RK** 
- [ ] Add **CataPro** report.  
- [ ] Consider **class-based** structure rather than independent functions ? 
- [ ] Move **PubChem API** queries to a dedicated module ? 

---
