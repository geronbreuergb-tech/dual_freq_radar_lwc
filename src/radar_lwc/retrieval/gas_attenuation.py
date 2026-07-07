"""
Gas attenuation coefficient α(f, h) [dB/km] from ITU-R P.676.

Uses the `itur` package (ITU-R P.676-13). Same physics as MATLAB's gaspl.
"""

import numpy as np
import xarray as xr
import itur
import astropy.units as u

def specific_gas_attenuation(freq_GHz, P_hPa, T_K, rho_v):
    """
    One-way specific gas attenuation from oxygen + water vapour.

    Parameters
    ----------
    freq_GHz : float
        Frequency in GHz (e.g. 35.0 or 94.0).
    P_hPa : float or array-like
        TOTAL atmospheric pressure in hPa. ITU expects dry+wet pressure.
    T_K : float or array-like
        Temperature in Kelvin.
    rho_v : float or array-like
        Water vapour density (absolute humidity) in g/m³.

    Returns
    -------
    gamma : same shape as inputs
        Specific attenuation in dB/km.
    """
    # itur returns astropy quantities; strip units with .value
    # Signature: gaseous_attenuation_inclined_path is for slant paths;
    # for specific attenuation at a point we use gammaw + gamma0.
   
    f  = freq_GHz * u.GHz
    P = P_hPa * u.hPa
    T = T_K * u.K
    rho = rho_v * u.g / u.m**3

    # Dry-air (oxygen) specific attenuation
    gamma_o = itur.models.itu676.gamma0_exact(f, P, rho, T).value
    # Water-vapour specific attenuation
    gamma_w = itur.models.itu676.gammaw_exact(f, P, rho, T).value
    return gamma_o + gamma_w


def alpha_field(
    P_field: xr.DataArray,
    T_field: xr.DataArray,
    AH_field: xr.DataArray,
    freq_ghz: float,
) -> xr.DataArray:
    """
    Vectorized gaseous attenuation α over a (Time, range) field.
    """

    # Convert xarray objects to NumPy arrays
    P_np = P_field.values
    T_np = T_field.values
    AH_np = AH_field.values

    # Compute gaseous attenuation
    alpha_np = specific_gas_attenuation(
        freq_ghz,
        P_np,
        T_np,
        AH_np,
    )

    return xr.DataArray(
        alpha_np,
        dims=T_field.dims,
        coords=T_field.coords,
        name=f"alpha_{freq_ghz:.0f}GHz",
        attrs={
            "units": "dB km-1",
            "long_name": f"Gaseous attenuation at {freq_ghz} GHz",
            "model": "ITU-R P.676",
        },
    )
