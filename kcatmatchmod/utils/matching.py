import logging
import pandas as pd 
import numpy as np

from kcatmatchmod.utils.temperature import arrhenius_equation
from kcatmatchmod.utils.organism import closest_org


# --- Check parameters ---


def check_enzyme(candidate, kcat_dict): 
    """
    Checks whether the enzyme in a candidate entry matches the model's enzyme.
    TODO: Manage to find the catalytic enzyme if multiple UniProt IDs are provided
    """
    if pd.notna(kcat_dict['uniprot_model']):
        enzyme_model = [e.strip() for e in kcat_dict.get('uniprot_model', '').split(';') if e.strip()]
        if candidate["UniProtKB_AC"] in enzyme_model:
            return 0
    return 2


def check_organism(candidate, general_criteria): 
    """
    Checks whether the organism in a candidate entry matches the expected organism.
    """
    if candidate["Organism"] == general_criteria["Organism"]:
        return 0
    return 2


def check_variant(candidate):
    """
    Checks whether the enzyme variant in a candidate entry is wildtype or unknown.
    """
    if candidate["EnzymeVariant"] == "wildtype":
        return 0
    else:  # Unknown
        return 1


def check_pH(candidate, general_criteria):
    """
    Checks whether the pH in a candidate entry matches the expected pH.
    """
    ph_min, ph_max = general_criteria["pH"]
    candidate_ph = candidate.get("pH", None)
    if candidate_ph is None:
        return 1 
    elif ph_min <= candidate_ph <= ph_max:    
        return 0
    else:  # Out of range 
        return 2
    

def check_temperature(candidate, general_criteria, api_output): 
    """
    Checks whether the temperature in a candidate entry matches the expected temperature.
    If the temperature is within the specified range is not met, verify if the Arrhenius equation can be applied.
    """
    temp_min, temp_max = general_criteria["Temperature"]
    candidate_temp = candidate.get("Temperature", None)
    if candidate_temp is not None and temp_min <= candidate_temp <= temp_max:
        return 0, False 
    ph_min, ph_max = general_criteria["pH"]
    filt = (
        api_output["pH"].between(ph_min, ph_max)
        & (api_output["UniProtKB_AC"] == candidate["UniProtKB_AC"])
        & api_output["Temperature"].notna()
    )
    temps_dispo = api_output.loc[filt, "Temperature"].nunique()
    if temps_dispo >= 2:
        # Go Arrhenius
        return 0, True
    elif candidate_temp is None:
        return 1, False
    else:
        return 2, False


def check_substrate(candidate, kcat_dict):
    """
    Checks whether the substrate in a candidate entry matches the model's substrates.
    """
    api = candidate['db']

    if api == 'sabio_rk':
        # 1. Try to match using KEGG reaction ID if provided
        if kcat_dict['KEGG_rxn_id'] == candidate['KeggReactionID']:
            # 1.1 Look for a substrate match
            model_substrates = [s.strip().lower() for s in kcat_dict.get('substrates_name', '').split(';') if s.strip()]
            candidate_substrates = [s.strip().lower() for s in candidate['Substrate'].split(';') if s.strip()]
            if any(substrate in candidate_substrates for substrate in model_substrates):
                return 0
            # 1.2 Look for a product match
            model_products = [p.strip().lower() for p in kcat_dict.get('products_name', '').split(';') if p.strip()]
            candidate_products = [p.strip().lower() for p in candidate['Product'].split(';') if p.strip()]
            if any(product in candidate_products for product in model_products):
                return 0
        # 2. Match using substrate names if no KEGG
        model_substrates = [s.strip().lower() for s in kcat_dict.get('substrates_name', '').split(';') if s.strip()]
        candidate_substrates = [s.strip().lower() for s in candidate['Substrate'].split(';') if s.strip()]
        if any(substrate in candidate_substrates for substrate in model_substrates):
            return 0  
        return 3
        
    elif api == 'brenda':
        model_substrates = [s.strip().lower() for s in kcat_dict.get('substrates_name', '').split(';') if s.strip()]
        candidate_substrates = [s.strip().lower() for s in candidate['Substrate'].split(';') if s.strip()]
        if any(substrate in candidate_substrates for substrate in model_substrates):
            return 0  
        return 3
    
    else:
        logging.error(f"Unknown API in 'db' column: {api}. Supported: 'sabio_rk', 'brenda'.")
        

# --- Scoring ---


def compute_score(kcat_dict, candidate, general_criteria, api_output, best_score):
    """
    Compute a score for the candidate based on the Kcat dictionary and general criteria.
    """
    score = 0
    # Check enzyme
    score += check_enzyme(candidate, kcat_dict)
    # Check organism
    score += check_organism(candidate, general_criteria)
    # Check variant
    score += check_variant(candidate) 
    # Check pH
    score += check_pH(candidate, general_criteria)
    # Check substrate 
    score += check_substrate(candidate, kcat_dict)
    if score > best_score: 
        return np.inf, False
    # Check temperature 
    score_temp, arrhenius = check_temperature(candidate, general_criteria, api_output) 
    score += score_temp
    return score, arrhenius


# --- Main --- 


def find_best_match(kcat_dict, api_output, general_criteria):
    """
    Finds the best matching enzyme entry from the provided API output based on: 
        - Kcat specific criteria: 
            * Substrate 
            * Enzyme 
             
       - General criteria : 
           * Organism
           * Temperature
           * pH

    This function filters out mutant enzyme variants, orders the remaining entries based on enzyme and organism similarity,
    and iteratively computes a score for each candidate to identify the best match. If a candidate requires an Arrhenius
    adjustment, the kcat value is recalculated accordingly.

    Parameters:
        kcat_dict (dict): Dictionary containing enzyme information.
        api_output (pd.DataFrame): DataFrame containing kcat entries and metadata from an API.
        general_criteria (dict): Dictionary specifying matching criteria.

    Returns:
        tuple:
            best_score (float): The lowest score found, representing the best match.
            best_candidate (dict or None): Dictionary of the best matching candidate's data, or None if no match is found.
    """
    # 1. Remove mutant
    api_output = api_output[api_output["EnzymeVariant"].isin(['wildtype', None])]

    # 2. Order based on the enzyme and organism
    api_output = closest_org(kcat_dict, api_output)

    # 3. Find the best match
    best_score = np.inf
    best_candidate = None

    for candidate in api_output.itertuples():
        candidate_dict = candidate._asdict()
        # Compute the score for the candidate
        score, arrhenius = compute_score(kcat_dict, candidate_dict, general_criteria, api_output, best_score)
        if score < best_score:
            if arrhenius: 
                kcat = arrhenius_equation(candidate_dict, api_output, general_criteria)
                candidate_dict['kcat'] = kcat 
                candidate_dict['Temperature'] = np.mean(general_criteria["Temperature"])  
            best_score = score
            best_candidate = candidate_dict
            if score == 0:  # Perfect match found
                break

    return best_score, best_candidate