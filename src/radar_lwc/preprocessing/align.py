import xarray as xr
import pandas as pd

def align_to_reference_time(
    da_ref: xr.DataArray,
    da_target: xr.DataArray,
    method: str = "nearest",
    tolerance: pd.Timedelta = pd.Timedelta("1s"),
) -> xr.DataArray:
    """
    Align a target DataArray to the time coordinate of a reference DataArray.

    Useful in radar intercomparisons where one instrument (e.g. Ka-band)
    is chosen as the time reference and the other (e.g. W-band) is
    matched to its sampling times. This preserves the reference radar's
    original timestamps and avoids creating artificial time grids.

    Parameters
    ----------
    da_ref : xr.DataArray
        Reference DataArray whose 'Time' coordinate defines the target grid.
    da_target : xr.DataArray
        DataArray to be re-sampled onto da_ref's time coordinate.
    method : {'nearest', 'linear'}, default 'nearest'
        - 'nearest' : nearest-neighbour matching (NO smoothing, preserves dBZ statistics).
                      Recommended for radar reflectivity comparisons.
        - 'linear'  : linear interpolation in time (smooths gradients slightly).
    tolerance : pd.Timedelta, default 1 second
        Maximum allowed time difference for nearest-neighbour matching.
        Samples without a partner within ± tolerance become NaN.
        Only applies if method='nearest'.

    Returns
    -------
    xr.DataArray
        da_target evaluated at da_ref['Time'].

    Notes
    -----
    For W-/Ka-band cloud radar comparisons sampling at ~1 Hz, a tolerance
    of 1 s is appropriate. Tighten it (e.g. 0.5 s) for stricter matching.
    """

    if method == "nearest":
        return da_target.reindex(
            Time=da_ref["Time"],
            method="nearest",
            tolerance=tolerance,
        )

    elif method == "linear":
        return da_target.interp(Time=da_ref["Time"])

    else:
        raise ValueError(f"method must be 'nearest' or 'linear', got {method!r}")
    
