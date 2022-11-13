"""
Utility functions for the grid
"""
import numpy as np

def json_encoder(obj):
    """
    To handle output of new types
    """
    if type(obj).__module__ == np.__name__:
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return obj.item()
    else:

        try:
            string_version = str(obj)
            return string_version
        except:
            raise TypeError(
                "Unserializable object {} of type {}. Attempted to convert to string but that failed.".format(
                    obj, type(obj)
                )
            )

    raise TypeError("Unknown type:", type(obj))
