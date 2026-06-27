import numpy as np
import xarray as xr

def alpha_gas_field_placeholder(template: xr.DataArray, freq_ghz: float) -> xr.DataArray:
    """
    PLACEHOLDER for one-way gas attenuation [dB/km].

    Replace with ITU-R P.676 / Rosenkranz / pyrtlib when
    pressure, humidity and temperature profiles are available.
    """
    if np.isclose(freq_ghz, 35.0):
        val = 0.10
    elif np.isclose(freq_ghz, 94.0):
        val = 1.00
    else:
        raise ValueError(f"No placeholder α implemented for {freq_ghz} GHz")

    return xr.DataArray(
        np.full(template.shape, val, dtype=float),
        dims=template.dims,
        coords=template.coords,
        name=f"alpha_{freq_ghz:.0f}GHz",
        attrs={
            "units": "dB km-1",
            "long_name": f"Gas attenuation at {freq_ghz} GHz (PLACEHOLDER)",
            "warning": "Constant placeholder; replace with real gas model.",
        },
    )