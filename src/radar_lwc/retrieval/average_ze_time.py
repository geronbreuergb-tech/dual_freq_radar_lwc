import numpy as np

def avg_ze_time(ze_db, window_s=30):
    """Average dBZ in linear units over a time window."""
    z_lin = 10.0 ** (ze_db / 10.0)          # Make Reflectivity Linear again

    # xarray rolling on the Time dim; center=True keeps the same time axis
    # Convert 'window_s' to a rolling count using the median dt
    dt = float(np.median(np.diff(ze_db["Time"].values).astype("timedelta64[s]").astype(float)))    #median time between "Time".values of the reflectivity data array in seconds (as float)
    n = max(1, int(round(window_s / dt)))         # Window size divided by average time step gives the number of time steps to include in the rolling average
    z_avg = z_lin.rolling(Time=n, center=True, min_periods=1).mean()      # Takes middle point( center=True) of the rolling window and computes the mean of the linear reflectivity values in that window. min_periods=1 means that if there is at least one valid value in the window, it will compute the mean; otherwise, it will return NaN.
                                               #min_periods=1: Compute an average as long as there is at least one valid value in the window. At the beginning and end of the time series, the average is then calculated from the available values only
                                                # min_periods=n: Compute an average only if the full window of n values is available. Otherwise, the result is NaN.
    return 10.0 * np.log10(z_avg)