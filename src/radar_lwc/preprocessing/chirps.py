import xarray as xr

def combine_chirps(ds_zen, var_prefix="ZE"):
    """
    Combine the three chirp variables (C1, C2, C3) of an RPG cloud radar
    into a single continuous profile along the range dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        The loaded RPG LV1 dataset (with C1Range, C2Range, C3Range dims).
    var_prefix : str
        The variable name without the chirp prefix.
        Examples: 'ZE' → combines C1ZE, C2ZE, C3ZE.
                  'MeanVel' → combines C1MeanVel, C2MeanVel, C3MeanVel.

    Returns
    -------
    combined : xarray.DataArray
        Variable on a unified (Time, range) grid.
    """
    # 1. Grab the three chirp arrays (DataArrays)
    c1 = ds_zen[f"C1{var_prefix}"]
    c2 = ds_zen[f"C2{var_prefix}"]
    c3 = ds_zen[f"C3{var_prefix}"]

    # 2. Rename their per-chirp range dim to a common name 'range'
    c1 = c1.rename({"C1Range": "range"})
    c2 = c2.rename({"C2Range": "range"})
    c3 = c3.rename({"C3Range": "range"})

    # 3. Assign each chirp its real range values (in metres) as the coordinate
    c1 = c1.assign_coords(range=ds_zen["C1Range"].values)
    c2 = c2.assign_coords(range=ds_zen["C2Range"].values)
    c3 = c3.assign_coords(range=ds_zen["C3Range"].values)

    # 4. Concatenate along the range axis
    combined = xr.concat([c1, c2, c3], dim="range")

    return combined
