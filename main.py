"""
Main script for the L2 mass loss routines by Wenbin Lu

TODO: Import functions
TODO: for each mass value pair gridpoint evolve an interpolation table
TODO: set up exploration ensemble to find the correct values to use
TODO: generate interpolation table with all the values for binary_c
TODO: make generate_header function for binary_c to use as the interpolation
"""

import os
import pickle
import pandas as pd
import numpy as np

grid_settings = {
    ###########
    # Grid settings to control the parameter pairs that we run wenbins functions for
    'log10M_accretor_array': np.arange(-1, 2, 10),
    'log10massratio_donor_accretor_array': np.arange(-1, 1, 10),
    ###########
    # Parameter configuration for wenbins functions
    # Separation
    'log10separation_min': -5,
    'log10separation_max': 1,
    'log10separation_N': 10,
    # Donor Mdot
    'log10Mdot_donor_min': -5,
    'log10Mdot_donor_max': -2,
    'log10Mdot_donor_N': 10,
    ###########
    # Control parameters for wenbins functions
    "savedir": os.path.join(
        os.environ["PROJECT_DATA_ROOT"],
        "thick_disk_calculations_wenbin",
        "TEST",
    ),
}
