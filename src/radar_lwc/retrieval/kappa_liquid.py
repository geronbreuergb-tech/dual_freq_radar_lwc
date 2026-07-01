import numpy as np
from radar_lwc.retrieval.dielectric import epsilon_water_liebe1989
import xarray as xr

C_LIGHT  = 2.998e8   # m/s
RHO_LIQ  = 1.0e6     # g/m^3 (i.e., 1000 kg/m^3 in g/m^3)

def kappa_liquid(freq_ghz: float, T_celsius: float) -> float:
    """
    One-way liquid water specific attenuation coefficient (dB km^-1 (g m^-3)^-1)
    in the Rayleigh / small-particle regime (Doviak & Zrnic 1993).

    κ = 4.343e3 * (6π/λ) / ρ_l * Im[ -(ε-1)/(ε+2) ]

    Parameters
    ----------
    freq_ghz : float
    T_celsius : float

    Returns
    -------
    float
        κ in dB·km⁻¹·(g·m⁻³)⁻¹
    """
    eps   = epsilon_water_liebe1989(freq_ghz, T_celsius)
    lam_m = C_LIGHT / (freq_ghz * 1e9)            # wavelength in m
    K     = (eps - 1.0) / (eps + 2.0)            # dielectric factor (note +2, not +1!)
    
    # Note: original Zhu equation as written has (ε+1) — that's a TYPO in the paper.
    # The correct Clausius-Mossotti form is (ε-1)/(ε+2). Verified against
    # Doviak & Zrnic 1993 eq. 8.55, Hogan et al. 2005.

    
    
    kappa = 4.343e3 * (6.0 * np.pi / lam_m) * np.imag(-K) / RHO_LIQ
    
    return kappa   # dB/km per g/m^3



def kappa_field(T_field: xr.DataArray, freq_ghz: float) -> xr.DataArray:
    """
    Vectorized κ over a (time, range) temperature field.
    """
    kappa_np = kappa_liquid(freq_ghz, T_field.values)

    return xr.DataArray(
        kappa_np,
        dims=T_field.dims,
        coords=T_field.coords,
        name=f"kappa_{freq_ghz:.0f}GHz",
        attrs={
            "units": "dB km-1 (g m-3)-1",
            "long_name": f"Rayleigh attenuation coefficient at {freq_ghz} GHz",
            "model": "Liebe et al. 1989 double-Debye",
        },
    )
