import logging 
from cobra.io import load_json_model, load_matlab_model, read_sbml_model


def read_model(model_path: str):
    """
    Reads a metabolic model from a given path.
    
    Parameters:
    - model_path: str
        Path to a model file.
    
    Returns:
    - model: COBRA.Model
        The COBRA model object.
    """
    # Check the file extension
    if model_path.endswith(".json"):
        return load_json_model(model_path)
    elif model_path.endswith(".mat"):
        return load_matlab_model(model_path)
    elif model_path.endswith(".xml") or model_path.endswith(".sbml"):
        return read_sbml_model(model_path)
    else:
        logging.error(f"Unsupported model file format: {model_path}")