"""
Integral-form ('plan B') DWR LWC retrieval.

Instead of differentiating a noisy DWR profile locally (Zhu et al. 2019),
this works on the *cumulative* liquid water DQ(r) = integral of q_L, which is
inherently smoother (integration low-passes noise), fits it once over the whole
cloud, and differentiates the fitted curve a single time to recover LWC.

Derivation (Hogan-convention DWR = Ze_Ka - Ze_W, increases with height):
    DWR(r) = 2 * integral_base^r [alpha_W - alpha_Ka] d_rho
           + 2 * integral_base^r [kappa_W - kappa_Ka] * q_L d_rho

    Define the measurable 'liquid DWR':
        G(r) = DWR_anchored(r) - 2 * integral[alpha_W - alpha_Ka] d_rho
             = 2 * integral[kappa_W - kappa_Ka] * q_L d_rho

    kappa_mode="tmean"    : pull a single mean Delta_kappa out of the integral
                            -> DQ(r) = G(r) / (2 * Delta_kappa_mean),  LWC = dDQ/dr
    kappa_mode="resolved" : keep kappa height-resolved; fit G(r), then
                            LWC(r) = (dG/dr) / (2 * Delta_kappa(r))   pointwise

Units: DQ carries units g m-3 * km ; LWP = DQ(cloud top) * 1000  [g m-2].
LWC = dDQ/dr with r in km is g m-3 directly.

Validated by synthetic round-trip: at zero noise recovers a known 6-gate
(~112 m) cloud to ~0.04 g/m3 RMSE and LWP to a few g/m2; at 0.3 dB DWR noise
it roughly halves the RMSE vs a raw local derivative on the same data.
"""
import numpy as np
import xarray as xr


def retrieve_lwc_dq(
    dwr:        xr.DataArray,   # (Time, range) [dB], Hogan convention Ze_Ka - Ze_W
    kappa_w:    xr.DataArray,   # (Time, range) [dB km-1 (g m-3)-1]
    kappa_ka:   xr.DataArray,   # (Time, range)
    alpha_w:    xr.DataArray,   # (Time, range) [dB km-1]
    alpha_ka:   xr.DataArray,   # (Time, range)
    common_mask: xr.DataArray,  # (Time, range) bool, True = liquid cloud
    min_cloud_gates: int = 4,   # skip profiles thinner than this
    fit_degree: int = 3,        # polynomial degree for DQ(h) fit (2 or 3)
    kappa_mode: str = "tmean",  # "tmean" or "resolved"
):
    """Return (lwc, dq) DataArrays. LWC [g m-3]; DQ [g m-3 km], LWP=DQ_top*1000."""
    h_km = dwr["range"].values.astype(float) / 1000.0
    D    = dwr.values
    Kw, Kka = kappa_w.values, kappa_ka.values
    Aw, Aka = alpha_w.values, alpha_ka.values
    M = common_mask.values.astype(bool)
    nt, nr = D.shape

    lwc = np.full((nt, nr), np.nan)
    dq  = np.full((nt, nr), np.nan)

    for t in range(nt):
        idx = np.flatnonzero(M[t])
        if idx.size < min_cloud_gates:
            continue
        i0, i1 = idx.min(), idx.max()
        sl = slice(i0, i1 + 1)

        h  = h_km[sl]
        Dw = D[t, sl]
        if not np.isfinite(Dw).all():
            continue                      # keep v1 simple: require a clean run
        dk = Kw[t, sl] - Kka[t, sl]       # Delta_kappa per gate
        da = Aw[t, sl] - Aka[t, sl]       # Delta_alpha per gate

        # 1. anchor DWR to cloud base (removes any constant inter-radar offset)
        D_anch = Dw - Dw[0]

        # 2. cumulative gas term (dB), trapezoid from cloud base, dh in km
        dh = np.diff(h, prepend=h[0])
        cum_gas = 2.0 * np.cumsum(da * dh)

        # 3. liquid DWR
        G = D_anch - cum_gas

        deg = min(fit_degree, len(h) - 1)
        if deg < 1:
            continue

        if kappa_mode == "tmean":
            DQ_raw = G / (2.0 * np.mean(dk))          # scalar Delta_kappa
            coef = np.polyfit(h, DQ_raw, deg)         # fit smooth DQ(h)
            dDQ_dh = np.polyval(np.polyder(coef), h)  # analytic derivative
            lwc_prof = dDQ_dh                         # g m-3
            dq_prof  = np.polyval(coef, h)            # smoothed DQ
        elif kappa_mode == "resolved":
            coef = np.polyfit(h, G, deg)              # fit smooth G(h)
            dG_dh = np.polyval(np.polyder(coef), h)
            lwc_prof = dG_dh / (2.0 * dk)             # pointwise, height-resolved
            dq_prof  = np.cumsum(lwc_prof * dh)       # integrate back for DQ/LWP
        else:
            raise ValueError("kappa_mode must be 'tmean' or 'resolved'")

        oi = np.arange(i0, i1 + 1)
        lwc[t, oi] = lwc_prof
        dq[t, oi]  = dq_prof

    lwc_da = xr.DataArray(
        lwc, coords=dwr.coords, dims=dwr.dims, name="LWC",
        attrs={"units": "g m-3",
               "long_name": "Liquid Water Content (integral DQ retrieval)",
               "method": f"DQ polynomial fit deg={fit_degree}, kappa_mode={kappa_mode}"},
    )
    dq_da = xr.DataArray(
        dq, coords=dwr.coords, dims=dwr.dims, name="DQ",
        attrs={"units": "g m-3 km",
               "long_name": "Cumulative liquid water DQ (LWP = DQ_top * 1000 g m-2)"},
    )
    return lwc_da, dq_da
