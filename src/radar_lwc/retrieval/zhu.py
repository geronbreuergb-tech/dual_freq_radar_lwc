"""
Zhu et al. (2019) sliding-window DWR-gradient LWC retrieval.

Reference
---------
Zhu, Z., Kollias, P., & Yang, F. (2019).
'The vertical structure of liquid water content in shallow clouds as
retrieved from dual-wavelength radar observations.'
JGR Atmospheres, 124, 14184–14197.
"""
import numpy as np
import xarray as xr


def retrieve_lwc_zhu(
    dwr:        xr.DataArray,       # (time, range)  [dB]
    kappa_w:    xr.DataArray,       # (time, range)  [dB·km⁻¹·(g·m⁻³)⁻¹]
    kappa_ka:   xr.DataArray,       # (time, range)
    alpha_w:    xr.DataArray,       # (time, range)  [dB·km⁻¹]
    alpha_ka:   xr.DataArray,       # (time, range)
    cloud_mask: xr.DataArray,       # (time, range)  bool — True = liquid cloud
    N_default:  int = 9,            # window size (gates)
    min_cloud_gates: int = 3,       # skip profiles with fewer than this
) -> xr.DataArray:
    """
    Retrieve LWC(time, range) via sliding-window polynomial DWR gradient.
    
    Algorithm (per time step):
    1. Identify cloud gates from cloud_mask.
    2. If cloud thickness < N_default, set N = thickness (Zhu's thin-cloud rule).
       If thickness < min_cloud_gates, skip profile entirely.
    3. Slide a window of N gates from cloud base to cloud top, step = 1 gate.
    4. In each window: fit DWR(h) ~ a·h² + b·h + c → dDWR/dh = 2ah + b.
    5. LWC[h] = (½·dDWR/dh − Δα(h)) / Δκ(h)   — per-gate, no averaging of κ.
    6. Each gate accumulates LWC contributions from every window it belongs to.
    7. Final LWC[h] = mean of all accumulated estimates.

    Parameters
    ----------
    dwr, kappa_w, kappa_ka, alpha_w, alpha_ka : xr.DataArray
        All have dims (time, range) and share the same range axis [meters].
    cloud_mask : xr.DataArray of bool
        True where liquid cloud is present.
    N_default : int
        Default window size in gates. Adjust per Δh to keep ~200 m physical window.
    min_cloud_gates : int
        Skip cloud profiles with fewer contiguous cloud gates than this.

    Returns
    -------
    xr.DataArray  LWC(time, range) [g/m³]. NaN outside cloud.
    """
    # ── Pull NumPy views for speed ───────────────────────────────────────────
    r_m     = dwr["range"].values.astype(float)         # heights [m]
    n_t, n_r = dwr.shape
    
    dwr_v   = dwr.values
    kw_v    = kappa_w.values
    kka_v   = kappa_ka.values
    aw_v    = alpha_w.values
    aka_v   = alpha_ka.values
    mask_v  = cloud_mask.values.astype(bool)
    
    dkappa = kw_v - kka_v       # (time, range)  per-gate Δκ
    dalpha = aw_v - aka_v       # (time, range)  per-gate Δα
    
    # ── Output accumulators ──────────────────────────────────────────────────
    lwc_sum = np.zeros((n_t, n_r), dtype=float)
    lwc_cnt = np.zeros((n_t, n_r), dtype=int)
    

    derivative_sum = np.zeros((n_t, n_r), dtype=float)
    derivative_count = np.zeros((n_t, n_r), dtype=int)


    # ── Outer loop: time step by time step ───────────────────────────────────
    for t in range(n_t):
        # ── Find contiguous cloud segments at this time ─────────────────────
        # For simplicity (v1): treat the union of all cloud gates as one segment.
        # If your data has multi-layer clouds, you'll want to split into
        # contiguous runs and process each separately — see comment below.
        cloud_idx = np.where(mask_v[t])[0]                        # Makes Array of all the indices where the cloud mask is True (i.e., where there is a cloud) ; [0] is for the first element of the tuple returned by np.where, which is the array of indices
        if cloud_idx.size < min_cloud_gates:                       
            continue
        
        i_base = cloud_idx.min()        # cloud-base gate index
        i_top  = cloud_idx.max()        # cloud-top gate index
        thickness = i_top - i_base + 1
        
        # ── Adapt window size for thin clouds ────────────────────────────────
        N = min(N_default, thickness)                      # Searches Minimum between the default window size and the thickness of the cloud. If the cloud is thinner than the default window size, it will use the thickness as the window size.
        if N < 3:                       # need ≥ 3 points for a quadratic fit
            continue
        
        # ── Slide window from cloud base to cloud top ────────────────────────
        for i_start in range(i_base, i_top - N + 2):                 # Python excludes the last index, so we need to add 2 to include the last window that starts at i_top - N + 1
            i_end = i_start + N         # exclusive
            
            h_win    = r_m[i_start:i_end]
            dwr_win  = dwr_v[t, i_start:i_end]
            mask_win = mask_v[t, i_start:i_end]
            
            # Require ALL gates in window to be valid cloud (Zhu-style)
            if not mask_win.all():
                continue
            if np.any(~np.isfinite(dwr_win)):                  # If one of the values in the dwr_win is not finite, skip this window; ~ is used to invert the boolean array returned by np.isfinite, so that it returns True for non-finite values and False for finite values. np.any then checks if any of the values in the inverted array are True, which means that there is at least one non-finite value in dwr_win.
                continue
            
            # ── Fit quadratic & take analytical derivative ───────────────────
            #   DWR(h) ≈ a·h² + b·h + c   →   dDWR/dh = 2a·h + b
            a, b, _ = np.polyfit(h_win, dwr_win, 2)            # h_win and dwr_win are the x and y values for the polynomial fit, respectively. The 2 indicates that we want to fit a polynomial of degree 2 (i.e., a quadratic). np.polyfit returns the coefficients of the polynomial in descending order, so a is the coefficient of h², b is the coefficient of h, and _ is the constant term (which we don't need for the derivative).
            dDWR_dh_per_m = 2.0 * a * h_win + b      # dB / m
            dDWR_dh       = dDWR_dh_per_m * 1000.0   # → dB / km   ★ unit fix
            
            # ── Per-gate Δκ, Δα for THIS window (no averaging!) ──────────────
            dk_win = dkappa[t, i_start:i_end]
            da_win = dalpha[t, i_start:i_end]
            
            # ── Zhu eq. (1) ─────────────────────────────────────────────────
            #   LWC = (½ · dDWR/dh  − (α_W − α_Ka)) / (κ_W − κ_Ka)
            #   Robust to Δκ ≈ 0 (shouldn't happen for liquid, but safe):
            with np.errstate(divide="ignore", invalid="ignore"):           # If division is invalid just ignore it and return NaN for that gate. This is important because if Δκ is very small, the division could result in a very large value or NaN, which would skew the results. By ignoring the error, we can safely compute LWC for the valid gates and leave the invalid ones as NaN.
                lwc_win = (0.5 * dDWR_dh - da_win) / dk_win
            
            # ── Accumulate into running sum / count ──────────────────────────
            valid = np.isfinite(lwc_win)                                #Looks if the values in the lwc_win array are finite (i.e., not NaN or Inf). Returns a boolean array of the same shape as lwc_win, where True indicates a finite value and False indicates a non-finite value.
            idx = np.arange(i_start, i_end)

            lwc_sum[t, idx[valid]] += lwc_win[valid]
            lwc_cnt[t, idx[valid]] += 1


            derivative_sum[t, idx[valid]] += dDWR_dh
            derivative_count[t, idx[valid]] += 1
    
    # ── Average across all windows that contributed ──────────────────────────
    with np.errstate(invalid="ignore", divide="ignore"):          # If division is invalid just ignore it and return NaN for that gate. 
        lwc = np.where(lwc_cnt > 0, lwc_sum / lwc_cnt, np.nan)

    with np.errstate(invalid="ignore", divide="ignore"):
        derivative = np.where(derivative_count > 0, derivative_sum / derivative_count, np.nan)


    # # ── Wrap into xarray ─────────────────────────────────────────────────────
    # return xr.DataArray(
    #     lwc,
    #     dims=dwr.dims,
    #     coords=dwr.coords,
    #     name="LWC",
    #     attrs={
    #         "units": "g m-3",
    #         "long_name": "Liquid Water Content (Zhu et al. 2019 retrieval)",
    #         "method": "Sliding-window quadratic DWR gradient",
    #         "window_size_default_gates": N_default,
    #         "min_cloud_gates": min_cloud_gates,
    #     }, 
    # ), derivative


    lwc_da = xr.DataArray(
        lwc,
        dims=dwr.dims,
        coords=dwr.coords,
        name="LWC",
        attrs={
        "units": "g m-3",
        "long_name": "Liquid Water Content (Zhu et al. 2019 retrieval)",
        "method": "Sliding-window quadratic DWR gradient",
        "window_size_default_gates": N_default,
        "min_cloud_gates": min_cloud_gates,
        },
    )

    derivative_da = xr.DataArray(
        derivative,
        dims=dwr.dims,
        coords=dwr.coords,
        name="dDWR_dh",
        attrs={
        "units": "dB m-1",
        "long_name": "Vertical derivative of DWR",
        },
    )

    return lwc_da, derivative_da