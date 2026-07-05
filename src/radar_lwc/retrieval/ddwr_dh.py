import numpy as np
import xarray as xr

def compute_ddwr_dh(
    dwr: xr.DataArray,              # (time, range) [dB]
    cloud_mask: xr.DataArray,       # (time, range) bool
    N_default: int = 9,             # window size (gates)
    min_cloud_gates: int = 3,       # skip profiles with fewer than this
) -> xr.DataArray:
    """
    Same Algorithm as Zhu
    """

    # ── Pull NumPy views for speed ───────────────────────────────────────────
    r_m = dwr["range"].values.astype(float)      # heights [m]
    n_t, n_r = dwr.shape

    dwr_v = dwr.values
    mask_v = cloud_mask.values.astype(bool)

    # ── Output accumulators ──────────────────────────────────────────────────
    derivative_sum = np.zeros((n_t, n_r), dtype=float)
    derivative_count = np.zeros((n_t, n_r), dtype=int)

    # ── Outer loop: time step by time step ───────────────────────────────────
    for t in range(n_t):

        # ── Find cloud gates at this time step ───────────────────────────────
        cloud_idx = np.where(mask_v[t])[0]       # indices where cloud_mask is True

        if cloud_idx.size < min_cloud_gates:
            continue

        i_base = cloud_idx.min()
        i_top = cloud_idx.max()
        thickness = i_top - i_base + 1

        # ── Adapt window size for thin clouds ────────────────────────────────
        N = min(N_default, thickness)

        if N < 3:
            continue

        # ── Slide window from cloud base to cloud top ────────────────────────
        for i_start in range(i_base, i_top - N + 2):

            i_end = i_start + N

            h_win = r_m[i_start:i_end]
            dwr_win = dwr_v[t, i_start:i_end]
            mask_win = mask_v[t, i_start:i_end]

            # Require all gates in the window to be cloud
            if not mask_win.all():
                continue

            # Skip windows containing NaN or Inf
            if np.any(~np.isfinite(dwr_win)):
                continue

            # ── Fit quadratic & compute analytical derivative ────────────────
            # DWR(h) ≈ a·h² + b·h + c
            # dDWR/dh = 2ah + b

            a, b, _ = np.polyfit(h_win, dwr_win, 2)

            dDWR_dh_per_m = 2.0 * a * h_win + b
            dDWR_dh = dDWR_dh_per_m * 1000.0      # dB/km

            # ── Accumulate derivative estimates ──────────────────────────────
            idx = np.arange(i_start, i_end)

            derivative_sum[t, idx] += dDWR_dh
            derivative_count[t, idx] += 1

    # ── Average contributions from overlapping windows ──────────────────────
    with np.errstate(invalid="ignore", divide="ignore"):
        derivative = np.where(
            derivative_count > 0,
            derivative_sum / derivative_count,
            np.nan,
        )

    # ── Wrap into xarray ─────────────────────────────────────────────────────
    derivative_da = xr.DataArray(
        derivative,
        dims=dwr.dims,
        coords=dwr.coords,
        name="dDWR_dh",
        attrs={
            "units": "dB km-1",
            "long_name": "Vertical derivative of DWR",
            "method": "Sliding-window quadratic fit",
            "window_size_default_gates": N_default,
            "min_cloud_gates": min_cloud_gates,
        },
    )

    return derivative_da

