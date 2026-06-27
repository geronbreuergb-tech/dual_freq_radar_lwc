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
    N_default:  int = 7,            # window size (gates)
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
    
    # ── Outer loop: time step by time step ───────────────────────────────────
    for t in range(n_t):
        # ── Find contiguous cloud segments at this time ─────────────────────
        # For simplicity (v1): treat the union of all cloud gates as one segment.
        # If your data has multi-layer clouds, you'll want to split into
        # contiguous runs and process each separately — see comment below.
        cloud_idx = np.where(mask_v[t])[0]
        if cloud_idx.size < min_cloud_gates:
            continue
        
        i_base = cloud_idx.min()        # cloud-base gate index
        i_top  = cloud_idx.max()        # cloud-top gate index
        thickness = i_top - i_base + 1
        
        # ── Adapt window size for thin clouds ────────────────────────────────
        N = min(N_default, thickness)
        if N < 3:                       # need ≥ 3 points for a quadratic fit
            continue
        
        # ── Slide window from cloud base to cloud top ────────────────────────
        for i_start in range(i_base, i_top - N + 2):
            i_end = i_start + N         # exclusive
            
            h_win    = r_m[i_start:i_end]
            dwr_win  = dwr_v[t, i_start:i_end]
            mask_win = mask_v[t, i_start:i_end]
            
            # Require ALL gates in window to be valid cloud (Zhu-style)
            if not mask_win.all():
                continue
            if np.any(~np.isfinite(dwr_win)):
                continue
            
            # ── Fit quadratic & take analytical derivative ───────────────────
            #   DWR(h) ≈ a·h² + b·h + c   →   dDWR/dh = 2a·h + b
            a, b, _ = np.polyfit(h_win, dwr_win, 2)
            dDWR_dh_per_m = 2.0 * a * h_win + b      # dB / m
            dDWR_dh       = dDWR_dh_per_m * 1000.0   # → dB / km   ★ unit fix
            
            # ── Per-gate Δκ, Δα for THIS window (no averaging!) ──────────────
            dk_win = dkappa[t, i_start:i_end]
            da_win = dalpha[t, i_start:i_end]
            
            # ── Zhu eq. (1) ─────────────────────────────────────────────────
            #   LWC = (½ · dDWR/dh  − (α_W − α_Ka)) / (κ_W − κ_Ka)
            #   Robust to Δκ ≈ 0 (shouldn't happen for liquid, but safe):
            with np.errstate(divide="ignore", invalid="ignore"):
                lwc_win = (0.5 * dDWR_dh - da_win) / dk_win
            
            # ── Accumulate into running sum / count ──────────────────────────
            valid = np.isfinite(lwc_win)
            idx = np.arange(i_start, i_end)
            lwc_sum[t, idx[valid]] += lwc_win[valid]
            lwc_cnt[t, idx[valid]] += 1
    
    # ── Average across all windows that contributed ──────────────────────────
    with np.errstate(invalid="ignore", divide="ignore"):
        lwc = np.where(lwc_cnt > 0, lwc_sum / lwc_cnt, np.nan)
    
    # ── Wrap into xarray ─────────────────────────────────────────────────────
    return xr.DataArray(
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
