"""
Main script for the L2 mass loss routines by Wenbin Lu

TODO: Import functions
TODO: for each mass value pair gridpoint evolve an interpolation table
TODO: set up exploration ensemble to find the correct values to use
TODO: generate interpolation table with all the values for binary_c
TODO: make generate_header function for binary_c to use as the interpolation
"""

from l2_massloss_interpolation.L2massloss_fork.settings import grid_settings
from l2_massloss_interpolation.L2massloss_fork.fL2_grid import run_fl2_grid_for_gridpoint


def run_grid_main(settings, run_grid, generate_interpolation_table, generate_data_header):
    """
    Main function to run the grid
    """

    # Run grid with multiprocessing
    if run_grid:

        run_fl2_grid_for_gridpoint(settings)


    # Generate interpolation table based on the results of the grid
    if generate_interpolation_table:
        pass

    # Generate data header based on the results of the grid
    if generate_data_header:
        pass



if __name__=="__main__":
    run_grid_main(settings=grid_settings, run_grid=True, generate_interpolation_table=False, generate_data_header=False)