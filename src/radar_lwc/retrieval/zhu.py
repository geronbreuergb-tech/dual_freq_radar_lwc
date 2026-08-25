"""
Zhu et al. (2019) sliding-window DWR-gradient LWC retrieval.

Reference
---------
Zhu, Z., Kollias, P., & Yang, F. (2019).
'The vertical structure of liquid water content in shallow clouds as
retrieved from dual-wavelength radar observations.'
JGR Atmospheres, 124, 14184–14197.


This version additionally returns, for every gate, the EFFECTIVE SAMPLE SIZE
M_eff (Zhu Eq. 8), computed from the empirical autocorrelation of the
individual sliding-window LWC estimates that contribute to that gate. This is
what Eq. (7) needs to turn the single-window uncertainty into the final
averaged uncertainty:
 
    sigma_final(gate) = sigma_single(gate) / sqrt(M_eff(gate))

"""
import numpy as np
import xarray as xr


def _autocorrelation(estimates, max_lag):
    """
    Sample autocorrelation at lags 1..max_lag.
 
    estimates : 1-D array of the M sliding-window LWC estimates at ONE gate,
                ordered by the window that produced them.
    max_lag   : largest lag to evaluate (must be <= M-2 so each corrcoef sees
                at least two points).
    """
    correlations = []
    for lag in range(1, max_lag + 1):    # I go through all the lags
        shifted_later = estimates[lag:]     #
        shifted_early = estimates[:-lag]    
        # np.corrcoef returns a 2x2 matrix [[1, r], [r, 1]]; [0][1] is r.
        r = np.corrcoef(shifted_later, shifted_early)[0][1]    # For lag=1 (If lwc in one window is large, is the lwc of the window one shift later still large?) I compare the lag array withput the last entry and the lag array without the first entry. The difference between the two arrays at each point is that of two following points in the lag. For these arrays I look at the correlation. 
        # Afterwards I do the same for lag=2 (If lwc in one window is large, is the lwc of the window two shifts later still large?) which takes away the last two entries of the lag array and the first two entries of the lag array. Thats how I see the correlation between estimate 1 and estimate 3 and estimate 4 and 2.
        # I do this for every lag to see the correlation between window 1 and 2, between 1,3 between 1,4 ... In the end I sum this up then for M_eff
        correlations.append(r)
    return np.array(correlations)
 
 
def effective_sample_size(estimates, min_estimates_for_autocorr=4):
    """
    Zhu Eq. (8): M_eff = M / (1 + 2 * sum_k rho(k)).
 
    estimates : list/array of the M window LWC estimates at one gate.
 
    Returns M_eff, clipped to [1, M]:
      - cannot be less than 1 (that would mean averaging INCREASES error)
      - cannot exceed M   (that would mean more independent samples than we have)
 
    If there are too few estimates to measure rho(k) reliably, we fall back to
    M_eff = M (i.e. treat them as independent). With very few overlapping
    windows this makes little difference, and it avoids trusting a correlation
    estimated from 2-3 points.
    """
    estimates = np.asarray(estimates, dtype=float)
    estimates = estimates[np.isfinite(estimates)]
    n_estimates = len(estimates)
 
    if n_estimates < 2:
        return float(max(n_estimates, 1))
    if n_estimates < min_estimates_for_autocorr:
        return float(n_estimates)          # too few for a stable rho(k)
    if np.std(estimates) == 0:
        return float(n_estimates)          # all identical -> corrcoef undefined
 
    # Use lags up to M-2 so the shortest slice still has >= 2 points.    
    max_lag = n_estimates - 2     # Correlation measures whether two variables vary together. With only one observation there is no variance. NumPy therefore returns NaN so I have to stop at lag 7 for M/N=9
    rho = _autocorrelation(estimates, max_lag)
    rho = rho[np.isfinite(rho)]            # drop any NaN lags defensively
 
    denominator = 1.0 + 2.0 * np.sum(rho)
    m_eff = n_estimates / denominator
    return float(np.clip(m_eff, 1.0, n_estimates))









def retrieve_lwc_zhu(
    dwr:        xr.DataArray,       # (time, range)  [dB]
    kappa_w:    xr.DataArray,       # (time, range)  [dB·km⁻¹·(g·m⁻³)⁻¹]
    kappa_ka:   xr.DataArray,       # (time, range)
    alpha_w:    xr.DataArray,       # (time, range)  [dB·km⁻¹]
    alpha_ka:   xr.DataArray,       # (time, range)
    cloud_mask: xr.DataArray,       # (time, range)  bool — True = liquid cloud
    N_default:  int = 9,            # window size (gates)
    min_cloud_gates: int = 3,       # skip profiles with fewer than this
    max_ddwr_dh=50.0,      # dB/km — physical ceiling on |dDWR/dh|
    max_window_std=3.0,    # dB    — reject noisy windows
    edge_pad=1,            # skip edge gates near cloud base/top
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
    n_used = np.zeros(n_t, dtype=int)
    
    derivative_sum = np.zeros((n_t, n_r), dtype=float)
    derivative_count = np.zeros((n_t, n_r), dtype=int)
    

    # NEW: keep every individual window estimate for each gate, so we can
    # measure their autocorrelation afterwards. A dict keyed by (t, gate) is
    # memory-light because only cloudy gates get entries.
    
    window_estimates = {}   # {(t, gate_index): [estimate_from_window_1, ...]}


    # ── Outer loop: time step by time step ───────────────────────────────────
    for t in range(n_t):
        # ── Find contiguous cloud segments at this time ─────────────────────
        # For simplicity (v1): treat the union of all cloud gates as one segment.
        # If your data has multi-layer clouds, you'll want to split into
        # contiguous runs and process each separately — see comment below.
        cloud_idx = np.where(mask_v[t])[0]                        # Makes Array of all the indices where the cloud mask is True (i.e., where there is a cloud) ; [0] is for the first element of the tuple returned by np.where, which is the array of indices. (There could be other arrays in tuple if not)
        if cloud_idx.size < min_cloud_gates:                       
            continue
        
        i_base = cloud_idx.min()        # cloud-base gate index
        i_top  = cloud_idx.max()        # cloud-top gate index

        i_base += edge_pad
        i_top  -= edge_pad #remove edges 


        thickness = i_top - i_base + 1
        
        # ── Adapt window size for thin clouds ────────────────────────────────
        N = min(N_default, thickness)                      # Searches Minimum between the default window size and the thickness of the cloud. If the cloud is thinner than the default window size, it will use the thickness as the window size.
        if N < 3:                       # need ≥ 3 points for a quadratic fit
            continue
        n_used[t] = N
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
            if dwr_win.std() > max_window_std:    # guarded retrieval: if the standard deviation of the DWR values in the window exceeds the maximum allowed value, skip this window. This is a noise filter to avoid fitting a polynomial to a very noisy set of DWR values, which could lead to unphysical results.
                continue

            # ── Fit quadratic & take analytical derivative ───────────────────
            #   DWR(h) ≈ a·h² + b·h + c   →   dDWR/dh = 2a·h + b
            a, b, _ = np.polyfit(h_win, dwr_win, 2)            # h_win and dwr_win are the x and y values for the polynomial fit, respectively. The 2 indicates that we want to fit a polynomial of degree 2 (i.e., a quadratic). np.polyfit returns the coefficients of the polynomial in descending order, so a is the coefficient of h², b is the coefficient of h, and _ is the constant term (which we don't need for the derivative).
            dDWR_dh_per_m = 2.0 * a * h_win + b      # dB / m
            dDWR_dh       = dDWR_dh_per_m * 1000.0   # → dB / km   ★ unit fix

            if np.any(np.abs(dDWR_dh) > max_ddwr_dh):   #guarded retrieval: if the absolute value of the derivative exceeds the maximum allowed value, skip this window. This is a physical ceiling on the vertical gradient of DWR, which helps to avoid unphysical results due to noise or other issues.
                continue
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

            # NEW: store each valid per-gate estimate under its gate key.
            for gate_index, estimate_value in zip(idx[valid], lwc_win[valid]):
                window_estimates.setdefault((t, gate_index), []).append(estimate_value)     #adds the lwc_win value to the list of estimates for that gate in the window_estimates dictionary. If the key (t, gate_index) does not exist yet, it creates a new list for that key and appends the estimate value to it.
    
    # ── Average across all windows that contributed ──────────────────────────
    with np.errstate(invalid="ignore", divide="ignore"):          # If division is invalid just ignore it and return NaN for that gate. 
        lwc = np.where(lwc_cnt > 0, lwc_sum / lwc_cnt, np.nan)

    with np.errstate(invalid="ignore", divide="ignore"):
        derivative = np.where(derivative_count > 0, derivative_sum / derivative_count, np.nan)


    # -- NEW: per-gate effective sample size and window count ---------------
    m_eff = np.full((n_t, n_r), np.nan, dtype=float)    # Full array of NaN values with shape (n_t, n_r) to store the effective sample size for each gate. This will be filled in later with the computed values.
    n_windows = np.zeros((n_t, n_r), dtype=float)  # Full array of zeros with shape (n_t, n_r) to store the number of sliding windows that contributed to each gate. This will be filled in later with the computed values.
 
    for (t, gate_index), estimates in window_estimates.items():
        n_windows[t, gate_index] = len(estimates)     # Puts the estimates for each gate
        m_eff[t, gate_index] = effective_sample_size(estimates)   # Calculates the effective sample size for each gate


    # # ── Wrap into xarray ─────────────────────────────────────────────────────

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

    m_eff_da = xr.DataArray(
        m_eff, dims=dwr.dims, coords=dwr.coords, name="M_eff",
        attrs={
            "long_name": "Effective sample size per gate (Zhu et al. 2019 Eq. 8)",
            "note": "empirical autocorrelation of the per-window LWC estimates; "
                    "clipped to [1, M]. Divide sigma_single by sqrt(M_eff) for "
                    "the final averaged uncertainty (Eq. 7).",
        },
    )
 
    n_windows_da = xr.DataArray(
        n_windows, dims=dwr.dims, coords=dwr.coords, name="n_windows",
        attrs={"long_name": "Number of sliding windows contributing to each gate (M)"},
    )
    # with the other DataArrays
    n_used_da = xr.DataArray(n_used, dims=["Time"], coords={"Time": dwr["Time"]},
                         name="N_used",
                         attrs={"long_name": "window size actually used, gates"})


    return lwc_da, derivative_da, m_eff_da, n_windows_da, n_used_da