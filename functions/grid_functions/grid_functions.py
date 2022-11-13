"""
File containing functions to handle the grid of calculations. The specific functions to calculate all the results for each grid-point are in the repo L2massloss_fork
"""

import os
import json
import copy

import pickle

import numpy as np







# """
# Script containing function to convert the interpolation table file to the structure for binary_c

# The read-out file contains the summarized results of the ballistic trajectories. instead of first saving the trajectory data we now just process that immediately
# """

# import os

# import pandas as pd
# from ballistic_integrator.functions.integrator.trajectory_integrator_class import (
#     EXIT_CODES,
# )
# from ballistic_integrator.functions.lagrange_points.lagrange_points_sepinski import (
#     A_formula,
# )

# from ballistic_integration_code.scripts.L1_stream.functions.functions import (
#     filter_df_on_exit_codes,
#     get_git_info_and_time,
# )
# from ballistic_integration_code.scripts.L1_stream.stage_2.stage_settings import (
#     stage_2_settings,
# )


# def generate_data_header_stage_2(
#     input_interpolation_textfile, output_interpolation_textfile
# ):
#     """
#     Function to write the interpolation data as a header array
#     """

#     #####################
#     # Config
#     length_decimals = 6
#     input_parameter_list = [
#         "log10normalized_thermal_velocity_dispersion",
#         "A_factor",
#         "massratio_accretor_donor",
#     ]
#     output_parameter_list = [
#         "rmin",
#         "rcirc",
#         "fraction_onto_donor",
#         "angular_momentum_multiplier_self_accretion",
#     ]

#     # Read out interpolation textfile
#     df = pd.read_csv(input_interpolation_textfile, sep="\s+", header=0)

#     # Add columns to df
#     df["A_factor"] = A_formula(f=df["synchronicity_factor"])

#     ############
#     # Prepare for write-out

#     # Filter the dataframe on those that have a correct exit code
#     filtered_df, _, _ = filter_df_on_exit_codes(
#         df,
#         allowed_exit_codes=stage_2_settings["allowed_exit_codes"],
#         fail_silently=False,
#     )

#     # Make sure the order is correct
#     filtered_df = filtered_df[input_parameter_list + output_parameter_list]

#     # make sure to sort the columns
#     filtered_df = filtered_df.sort_values(by=input_parameter_list)

#     # Get the number of lines we will output to the table and some others
#     num_lines = len(filtered_df.index)
#     num_input_parameters = len(input_parameter_list)
#     num_output_parameters = len(output_parameter_list)

#     # Get indices of input parameters
#     index_log10normalized_soundspeed = input_parameter_list.index(
#         "log10normalized_thermal_velocity_dispersion"
#     )
#     index_A_factor = input_parameter_list.index("A_factor")
#     index_mass_ratio = input_parameter_list.index("massratio_accretor_donor")

#     # Get indices of output parameters
#     index_rmin = output_parameter_list.index("rmin")
#     index_rcirc = output_parameter_list.index("rcirc")
#     index_fraction_accretion_onto_donor = output_parameter_list.index(
#         "fraction_onto_donor"
#     )
#     index_angular_momentum_multiplier_self_accretion = output_parameter_list.index(
#         "angular_momentum_multiplier_self_accretion"
#     )

#     # Get git info
#     git_info = get_git_info_and_time()

#     with open(output_interpolation_textfile, "w") as f:
#         ########################
#         # Write top-level explanation
#         f.write(
#             "// Ballistic stream trajectory integration stage 2 dataset. This dataset aims to provide an interpolation dataset to determine rmin and racc for the stream based on the mass ratio q and the A factor that captures synchronicity, eccentricty and mean anomaly (see sepinksy 2007), to improve RLOF calculations and estimates for interactions and specifically for accretion disks.\n\n"
#         )

#         ########################
#         # Write git revision information
#         f.write(
#             "// Generated on: {} with git repository: {} branch: {} commit: {}\n\n".format(
#                 git_info["datetime_string"],
#                 git_info["repo_name"],
#                 git_info["branch_name"],
#                 git_info["commit_sha"],
#             )
#         )

#         ########################
#         # Write defines
#         f.write("#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_STAGE_2\n\n")
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_NUM_DATALINES {}\n\n".format(
#                 num_lines
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_NUM_INPUT_PARAMETERS {}\n\n".format(
#                 num_input_parameters
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_NUM_OUTPUT_PARAMETERS {}\n\n".format(
#                 num_output_parameters
#             )
#         )

#         ########################
#         # Add lines containing the input indices
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_INPUT_INDEX_LOG10NORMALIZED_SOUNDSPEED {}\n\n".format(
#                 index_log10normalized_soundspeed
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_INPUT_INDEX_A_factor {}\n\n".format(
#                 index_A_factor
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_INPUT_INDEX_MASS_RATIO {}\n\n".format(
#                 index_mass_ratio
#             )
#         )

#         ########################
#         # Add lines containing the output indices
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_OUTPUT_INDEX_RMIN {}\n\n".format(
#                 index_rmin
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_OUTPUT_INDEX_RCIRC {}\n\n".format(
#                 index_rcirc
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_OUTPUT_INDEX_FRACTION_SELF_ACCRETION {}\n\n".format(
#                 index_fraction_accretion_onto_donor
#             )
#         )
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_OUTPUT_INDEX_SPECIFIC_ANGULAR_MOMENTUM_CHANGE_FACTOR_SELF_ACCRETION {}\n\n".format(
#                 index_angular_momentum_multiplier_self_accretion  # TODO: this needs fixing!
#             )
#         )

#         # Add start for the data array
#         f.write(
#             "#define RLOF_HENDRIKS2022_BALLISTIC_STREAM_INTERPOLATION_TABLE_DATA \\\n"
#         )

#         # Write data to file
#         for lineno, line in enumerate(filtered_df.values.tolist()):
#             termination = ",\\\n" if not lineno == num_lines - 1 else "\n"
#             outputline = (
#                 ",".join([str(round(el, length_decimals)) for el in line]) + termination
#             )

#             f.write(outputline)


# if __name__ == "__main__":
#     # Get input file
#     input_filename = os.path.join(
#         os.getenv("PROJECT_DATA_ROOT"),
#         "ballistic_data/L1_stream/ballistic_stream_integration_results_stage_2",
#         "trajectory_summary_data.txt",
#     )

#     # set output file
#     if os.getenv("BINARY_C"):
#         output_dir = os.path.join(os.getenv("BINARY_C"), "src", "RLOF")

#     # output_dir = os.path.join(
#     #     os.getenv("PROJECT_DATA_ROOT"),
#     #     "ballistic_data/L1_stream/ballistic_stream_integration_results_stage_2",
#     # )

#     # output_dir = "results/"

#     #
#     output_filename = os.path.join(
#         output_dir, "RLOF_Hendriks2022_ballistic_stream_table.h"
#     )
#     os.makedirs(output_dir, exist_ok=True)

#     generate_data_header_stage_2(
#         input_interpolation_textfile=input_filename,
#         output_interpolation_textfile=output_filename,
#     )
