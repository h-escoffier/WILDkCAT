import logging
import pandas as pd 
import numpy as np


def arrhenius_equation(candidate, api_output, general_criteria):
    """
    Estimates the kcat value at a target temperature using the Arrhenius equation, based on available experimental data.

    Parameters:
        candidate (dict): Information about the enzyme candidate.
        api_output (pd.DataFrame): DataFrame containing experimental kcat values.
        general_criteria (dict): Dictionary specifying selection criteria, including 'Temperature' and 'pH'.

    Returns:
        float: Estimated kcat value at the objective temperature, calculated using the Arrhenius equation.
    """

    def find_ea(df, expected_range=(50000, 150000)):
        """
        Estimate the activation energy (Ea) using the Arrhenius equation from kcat values at different temperatures.

        Parameters:
            df (pd.DataFrame): DataFrame with at least 'Temperature' (°C) and 'value' (kcat) columns.
            expected_range (tuple): Expected range for Ea in J/mol. Default is (50000, 150000).

        Returns:
            float: Estimated activation energy (Ea) in J/mol. 
        """

        r = 8.314  # Gas constant in J/(mol*K)

        # Filter out rows with missing values
        valid = df[['Temperature', 'value']].dropna()

        # Convert temperature to Kelvin
        temps = valid['Temperature'].values
        temps_K = temps + 273.15
        kcats = pd.to_numeric(valid['value'], errors='coerce').values

        x = 1 / temps_K
        y = np.log(kcats)
        slope, _ = np.polyfit(x, y, 1)
        ea = float(-slope * r)

        if not (expected_range[0] <= ea <= expected_range[1]):
            logging.warning(f"Estimated Ea ({ea:.0f} J/mol) is outside the expected range {expected_range} J/mol.")

        return ea

    def calculate_kcat(temp_obj, ea, kcat_ref, temp_ref): 
        """
        Calculates the catalytic rate constant (kcat) at a given temperature using the Arrhenius equation.

        Parameters: 
            temp_obj (float): The target temperature (in Kelvin) at which to calculate kcat.
            ea (float): The activation energy calculated using find_ea(). 
            kcat_ref (float): The reference kcat value measured at temp_ref.
            temp_ref (float): The reference temperature (in Kelvin) at which kcat_ref was measured.

        Returns: 
            float: The calculated kcat value at temp_obj.
        """
        r = 8.314
        kcat_obj = kcat_ref * np.exp(ea / r * (1/temp_ref - 1/temp_obj))
        return kcat_obj

    # Objective temperature
    obj_temp = np.mean(general_criteria["Temperature"]) + 273.15

    # Format the api_output DataFrame
    ph_min, ph_max = general_criteria["pH"]
    filters = (
        (api_output["UniProtKB_AC"] == candidate["UniProtKB_AC"]) &
        api_output["Temperature"].notna() &
        api_output["value"].notna() &
        api_output["pH"].between(ph_min, ph_max)
    )
    api_filtered = api_output.loc[filters, ["Temperature", "value"]].copy()
    
    # Convert temperatures to Kelvin
    api_filtered["Temperature"] = api_filtered["Temperature"] + 273.15

    # Estimate the activation energy (Ea)
    ea = find_ea(api_filtered)

    # Select one kcat for the ref
    kcat_ref = float(api_filtered['value'].iloc[0])
    temp_ref = float(api_filtered['Temperature'].iloc[0])
        
    kcat = calculate_kcat(obj_temp, ea, kcat_ref, temp_ref)
    return kcat
