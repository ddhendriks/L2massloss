"""
Function containing routine to write the settings to a file
"""

import os
import json

from l2_massloss_interpolation.L2massloss_fork.functions.grid_functions.utils import json_encoder

def write_settings_to_file(settings, output_filename):
    """
    Function to write the settings to an output file
    """

    # Make dir if not existant
    dirname = os.path.dirname(output_filename)
    os.makedirs(dirname, exist_ok=True)

    # Write settings json
    with open(output_filename, "w") as f:
        f.write(json.dumps(settings, default=json_encoder))
