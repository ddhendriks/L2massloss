"""
Function containing the generator description
"""

import copy

def grid_generator(settings):
    """
    Generator function to generate gridpoint for the thickness calculations
    """

    # Set the parameter ranges
    log10M_accretor_array = settings["log10M_accretor_array"]
    log10massratio_donor_accretor_array = settings["log10massratio_donor_accretor_array"]

    # Loop over parameters and yield the settings
    for log10M_accretor in log10M_accretor_array:
        for log10massratio_donor_accretor in log10massratio_donor_accretor_array:
            # Copy settings
            settings_copy = copy.deepcopy(settings)

            # Write gridpoints: 
            settings_copy['mass_accretor'] = 10**log10M_accretor
            settings_copy['massratio_accretor_donor'] = 10**log10massratio_donor_accretor

            # Copy the settings
            returned_settings = copy.deepcopy(settings)

            yield returned_settings
