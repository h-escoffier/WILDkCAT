import pandas as pd 


def matching_function(kcat_dict, df, general_criterias, exclude_criteria_keys=None):
    pass 


def create_kcat_value(df, method='mean'):
    """
    Aggregate kcat values from the DataFrame using the specified method.

    Parameters:
        df (pd.DataFrame): DataFrame containing SABIO-RK entries.
        method (str): Aggregation method: 'mean', 'max', or 'min'.

    Returns:
        float or None: Aggregated kcat value, or None if no valid values.
    """
    # Ensure the column exists and is numeric
    if 'parameter.startValue' not in df.columns:
        return None
    values = pd.to_numeric(df['parameter.startValue'], errors='coerce').dropna()
    if values.empty:
        return None
    if method == 'mean':
        return values.mean()
    elif method == 'max':
        return values.max()
    elif method == 'min':
        return values.min()
    else:
        raise ValueError("Invalid 'method' parameter. Choose from 'mean', 'max', or 'min'.")