import xarray as xr
import numpy as np

def drop_duplicate_times(ds: xr.Dataset, time_dim: str = "Time") -> xr.Dataset:
    """
    Remove duplicate timestamps along the time dimension, keeping the first occurrence.
    Also sorts the dataset chronologically.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset that may contain duplicated time coordinates.
    time_dim : str
        Name of the time dimension. Default: 'Time'.

    Returns
    -------
    xarray.Dataset
        Cleaned dataset with unique, sorted timestamps.
    """
    # 1. Sort chronologically (safety)
    ds = ds.sortby(time_dim)
    
    # 2. Boolean mask: keep first occurrence of each timestamp
    _, unique_idx = np.unique(ds[time_dim].values, return_index=True)
    ds = ds.isel({time_dim: sorted(unique_idx)})
    
    return ds
