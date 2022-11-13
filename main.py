"""
Main script for the L2 mass loss routines by Wenbin Lu

TODO: set up exploration ensemble to find the correct values to use
TODO: generate interpolation table with all the values for binary_c
TODO: make generate_header function for binary_c to use as the interpolation
"""

import os

from l2_massloss_interpolation.L2massloss_fork.settings import grid_settings
from l2_massloss_interpolation.L2massloss_fork.functions.grid_functions.run_all import run_all


if __name__=="__main__":
    run_all(settings=grid_settings, run_grid=True, generate_interpolation_table=False, generate_data_header=False)