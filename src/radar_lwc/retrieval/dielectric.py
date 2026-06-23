import numpy as np

def epsilon_water_liebe1989(freq_ghz: float, T_celsius: float) -> complex:
    """
 
    
    Reference
    ---------
  Millimeter-wave attenuation and delay rates due to fog/cloud conditions    Liebe, Manabe , and Hufford, 1989, [In Zotero]
    Seite 2 von 7

    Parameters
    ----------
    freq_ghz : float
        Frequency in GHz.
    T_celsius : float
        Temperature in °C.

    Returns
    -------
    complex
        ε = ε' − j·ε''  

    Notes
    -----
    Liebe (1989) uses the convention

        ε = ε' - j ε''

    whereas many radar texts use

        ε = ε' + j ε''

    For attenuation calculations only the magnitude
    of the imaginary component is used.
    """
    # Inverse temperature variable (Liebe's convention)
    theta = 300.0 / (T_celsius + 273.15)
    tm1   = theta - 1.0                                 # (θ − 1)

    # Permittivity constants  — eq. (8c)
    eps0 = 77.66 + 103.3 * tm1
    eps1 = 5.48
    eps2 = 3.51

    # Relaxation frequencies (GHz) — eq. (8d)
    f_p  = 20.09 - 142.4 * tm1 + 294.0 * tm1**2
    f_s  = 590.0 - 1500.0 * tm1

    # Double-Debye complex permittivity — eq. (7)
    eps = (eps0 - eps1) / (1.0 + 1 * (freq_ghz / f_p)) \
        + (eps1 - eps2) / (1.0 + 1j * (freq_ghz / f_s)) \
        + eps2

    return eps   # ε' − j·ε''  (imaginary part is negative)
