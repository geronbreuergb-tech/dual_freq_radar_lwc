
import numpy as np

from radar_lwc.preprocessing.chirps import combine_chirps

def get_reflectivity(ds):
    """
    Combine the three chirps' linear reflectivity and convert to dBZ.

    Parameters
    ----------
    ds : xarray.Dataset
        Loaded RPG LV1 dataset.

    Returns
    -------
    xarray.DataArray
        Reflectivity Ze in dBZ on a (Time, range) grid.
    """
    # Combine the chirps (still in linear units mm^6/m^3)
    ze_lin = combine_chirps(ds, "ZE")
    
    # Mask fill value -999 and any non-positive values (cannot take log)
    ze_lin = ze_lin.where(ze_lin > 0)
    
    # Linear → dBZ
    ze_dbz = 10 * np.log10(ze_lin)
    
    # Add metadata
    ze_dbz.name = "Ze"
    ze_dbz.attrs["units"] = "dBZ"
    ze_dbz.attrs["long_name"] = "Equivalent reflectivity factor"
    
    return ze_dbz
