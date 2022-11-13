"""
Script containing the settings for the project
"""

import os
import numpy as np

this_file = os.path.abspath(__file__)
this_file_dir = os.path.dirname(this_file)

grid_settings = {
    ###########
    # Grid settings to control the parameter pairs that we run wenbins functions for
    'log10M_accretor_array': np.arange(-1, 2, 10),
    'log10massratio_donor_accretor_array': np.arange(-1, 1, 10),
    ###########
    # Parameter configuration for wenbins functions
    # Separation
    'log10separation_min': 0.15,
    'log10separation_max': 3.25,
    # logamin, logamax = 0.15, 3.25   # [Msun/yr]
    'log10separation_N': 10,
    # Na = 200

    # Donor Mdot
    'log10Mdot_donor_min': -5,
    'log10Mdot_donor_max': -1.3,
    # logM1dotmin, logM1dotmax = -5.3, -1.3   # [Rsun]
    'log10Mdot_donor_N': 10,
    # NM1dot = 200

    # Disk thickness search
    "disk_thickness_search_grid_min": 0.1,
    "disk_thickness_search_grid_max": 1.0,
    # thegrid_min, thegrid_max = 0.1, 1.   # 0.1 to 1 is sufficient, we use analytic result below 0.1
    "disk_thickness_search_grid_N": 10,
    # Nthe = 50   # ~50 is accurate enough

    ###########
    # Control parameters for wenbins functions
    "savedir": os.path.join(
        os.environ["PROJECT_DATA_ROOT"],
        "thick_disk_calculations_wenbin",
        "TEST",
    ),
    # savedir = './disk_sltns/'
    "nofL2": False,
    # nofL2 = False   # turn on/off fL2
    "kappa_case": 1, # 1 -- Hrich, solarZ; 2 -- Hpoor, solarZ; 3 -- Hrich, lowZ
    # case = 1
    "kappa_fdir": os.path.join(this_file_dir, 'data', 'kap_data'),
    # fdir = './kap_data/'
    "kappa_extrap": False,
    # extrap = False   # extrapolation beyond grid may not be accurate [no need to]
    "mass_accretor": 10.,
    #M2_in_Msun = 10.   # [Msun]
    "massratio_accretor_donor": .5,
    # q = .5     # mass ratio accretor/donnor = M2/M1
    "alpha_ss": 0.1,
    # alpha_ss = 0.1   # viscosity
    "tol": 1e-8,
    # tol = 1e-8   # fractional tolerance for bisection method
    "eps_small": 1e-12,
    # eps_small = 1e-12   # a very small number
    "Tfloor": 3e3,
    # Tfloor = 3e3   # [K] --- the minimum value for disk temperature solution
}
