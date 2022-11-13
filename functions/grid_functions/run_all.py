"""
Main function to run the grid
"""

import os

from david_phd_functions.backup_functions.functions import backup_if_exists

from l2_massloss_interpolation.L2massloss_fork.fL2_grid import run_fl2_grid_for_gridpoint

from l2_massloss_interpolation.L2massloss_fork.functions.grid_functions.multiprocessing_routines import multiprocessing_routine
from l2_massloss_interpolation.L2massloss_fork.functions.grid_functions.write_settings_to_file import write_settings_to_file
from l2_massloss_interpolation.L2massloss_fork.functions.grid_functions.grid_generator import grid_generator

def run_all(settings, run_grid, generate_interpolation_table, generate_data_header):
    """
    Main function to run the grid
    """


    ##############################################
    # Run simulations and handle processing of data
    if run_grid:

        ##############################################
        # Check if the directory of the current simulation exists. If it does, then delete the old dir
        savedir = settings["savedir"]
        if settings.get("backup_if_data_exists", True):
            backup_if_exists(savedir, remove_old_directory_after_backup=True)

        #####################
        # Write settings to file
        write_settings_to_file(
            settings,
            output_filename=os.path.join(settings["savedir"], "settings.json"),
        )

        ######################
        # Set up grid generator
        grid_generator_instance = grid_generator(settings)

        ######################
        # Call the multiprocessing routine and let it handle everything
        multiprocessing_routine(
            configuration=settings,
            job_iterable=grid_generator_instance,
            target_function=run_fl2_grid_for_gridpoint,
        )

    # Generate interpolation table based on the results of the grid
    if generate_interpolation_table:
        pass

    # Generate data header based on the results of the grid
    if generate_data_header:
        pass
