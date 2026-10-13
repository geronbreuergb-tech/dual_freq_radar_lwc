# VARIABLE AND NAME REFERENCE — dual-frequency radar LWC retrieval

Extracted directly from the notebooks, not from memory. Every name below was
read out of the `.ipynb` files or out of the RPG LV1 variable listing.

Companion to `PROJECT_HANDOFF.md`. That file explains the science; this file
lists the names. Paste this one into a new chat whenever the work touches
code, and paste the handoff when the work touches physics.

Source notebooks read for this file:
`02_align_w_ka`, `02_05_gas_attenuation`, `03_dwr_and_masking`,
`04_lwc_retrieval_zhu` (+ 4 variants), `04_zhu_guarded_retrieval_CLEANED`,
`05_lwc_retrieval_dq_integral`, `06_error_estimation`.

---

# PART A — FILE I/O MAP (exact paths, exact spellings)

`DS = Processed_Dataset = "20230930"` in every notebook.
Base path is always relative: `../data/processed/{DS}/`.

| NB | Reads | Writes |
|---|---|---|
| 02 | `../data/raw/cmtrace/{Dataset_94GHz}/*ZEN.LV1.nc`, `../data/raw/cmtrace/{Dataset_35GHz}/*ZEN.LV1.nc` | `{DS}/aligned_LWP2.nc` |
| 02_05 | `{DS}/aligned_LWP2.nc` | `{DS}/gas_attenuation.nc` |
| 03 | `{DS}/aligned_LWP2.nc`, `{DS}/gas_attenuation.nc` | `{DS}/retrieval_input.nc`, `{DS}/retrieval_input_4km.nc`, `{DS}/retrieval_input_4km_alpha_kappa_calculated.nc` |
| 04 CLEANED | `{DS}/retrieval_input_4km_alpha_kappa_calculated.nc` | `{DS}/retrieval_output_zhu.nc`, **`../data/processed/retrieval_output_zhu_pad0.nc`**, **`..._pad1.nc`** |
| 05 | `{DS}/retrieval_input_4km_alpha_kappa_calculated.nc` | `{DS}/retrieval_output_dq.nc` |
| 06 | `{DS}/retrieval_input_4km_alpha_kappa_calculated.nc`, `{DS}/retrieval_output_zhu.nc`, `{DS}/retrieval_output_dq.nc`, `../data/processed/retrieval_output_zhu_pad{0,1}.nc` | **`../data/processed/lwc_uncertainty.nc`**, **`../data/processed/lwc_uncertainty_pad{0,1}.nc`** |

**Bold paths have no `{DS}` subfolder.** They sit in `../data/processed/`
directly. Processing a second dataset overwrites them silently. Fix before
running PPI, or zenith and PPI results will be indistinguishable on disk.

### Dataset names

```python
Dataset_94GHz     = "20230930_94GHz"     # W-band raw folder
Dataset_35GHz     = "20230930_35GHz"     # Ka-band raw folder
Processed_Dataset = "20230930"
# the 28 Sept folders use a different convention:
#   "20230928"  (W)  and  "20230928_KaBand"  (Ka)
```

### Dead / superseded filenames (do not use, they appear in older notebooks)

```
retrieval_input_2km.nc                          # 04 variants; nb03 has a cell that DELETES it
retrieval_input_2km_alpha_kappa_calculated.nc   # 04_lwc_retrieval_zhu_common_alpha_kappa only
```

`retrieval_input_2km_alpha_kappa_calculated_scan.nc` is what the notebook
08/09 drafts open. **That name does not exist.** The current name is
`_4km_`. Either change the drafts or change notebook 03's output name.

### Naming trap: `_4km_` files are sliced to 4000 m, variables say `2km`

Notebook 03 slices with `range=slice(0, 4000)` but names the Python
variables `dwr_2km`, `cloud_mask_2km`, `ds_retrieval_2km`, etc. The `2km`
in the variable names is historical. The number that matters is **4000**.

---

# PART B — RPG LV1 RAW VARIABLES (complete, from the instrument listing)

`Time` is the last dimension in the file listing (`C1Range,Time`); xarray
exposes it as `(Time, C1Range)` after loading. Sizes below are from the
listing example: 3 chirps, 22/52/265 gates, 68 profiles.

### Per-chirp fields — prefix `C1`, `C2`, `C3`

| Name | Meaning | Units | Loaded now? |
|---|---|---|---|
| `C{n}ZE` | equivalent reflectivity factor | mm^6/m^3 linear | **yes** |
| `C{n}Range` | gate centres | m | **yes** |
| `C{n}MeanVel` | mean Doppler velocity | m/s | **yes** |
| `C{n}SLv` | sensitivity limit, vertical pol | mm^6/m^3 | **yes** |
| `C{n}SLh` | sensitivity limit, horizontal pol | mm^6/m^3 | **yes** |
| `C{n}SLDR` | **slanted linear depolarisation ratio** | dB | no — needed for the insect filter |
| `C{n}SpecWidth` | Doppler spectral width | m/s | no — replaces ASSUMPTION 2 in nb06 |
| `C{n}ZE45` | reflectivity, 45 deg slanted channel | mm^6/m^3 | no |
| `C{n}ZDR` | differential reflectivity | dB | no |
| `C{n}Skew`, `C{n}Kurt` | spectral skewness, kurtosis | - | no |
| `C{n}RHV` | co-cross channel correlation | - | no |
| `C{n}PhiDP`, `C{n}KDP` | differential phase, specific diff. phase | rad, deg/km | no |
| `C{n}DiffAtt` | differential attenuation | dB/km | no |
| `C{n}SRCX` | slanted correlation coefficient | - | no |

`_FillValue = -999` on the polarimetric fields.

### Per-profile, dims `(Time,)`

| Name | Meaning | Units | Loaded now? |
|---|---|---|---|
| `Time` | seconds since epoch 2001-01-01 | s | **yes** |
| `Timems` | millisecond part of `Time` | ms | no |
| **`Elv`** | **elevation** | **deg** | **no — REQUIRED for scanning** |
| **`Azm`** | **azimuth** | **deg** | **no — REQUIRED for scanning** |
| `TProf` | temperature profile (MWR) | K | **yes** |
| `AHProf` | absolute humidity profile | g/m3 | **yes** |
| `RHProf` | relative humidity profile | % | no |
| `LWP` | MWR liquid water path | g/m2 | **yes** |
| `SurfPres` | surface pressure | hPa | **yes** |
| `SurfTemp`, `SurfRelHum` | surface T, RH | K, % | no |
| `SurfWS`, `SurfWD` | surface wind speed, direction | m/s, deg | `SurfWS` used ad hoc in nb04 |
| `Rain` | rain flag / rate | - | no |
| `QualFlag` | quality flag | - | no |
| `Status` | heater / blower mitigation flags | - | no |
| `TPow`, `TTemp`, `RTemp`, `PCTemp` | transmitted power, transmitter / receiver / PC temperature | W, K | no |
| `PowIF` | IF power at ADC | uW | no |
| `DDVolt`, `DDTb` | direct-detection voltage, brightness temperature | V, K | no |
| `Res_1` | active/passive channel noise power ratio | - | no |
| `Res_2`, `Res_3` | reserved | - | no |

**`ElvVel` and `AzmVel` do not exist in this LV1 version.** The handoff
suggests checking for them; they are not in the file. Scan rate has to be
derived as `diff(Azm) / diff(Time)`.

### Altitude grids — dims `(TAlt,)` and `(HAlt,)`

| Name | Meaning | Units |
|---|---|---|
| `TAlts` | altitudes of `TProf` | m |
| `HAlts` | altitudes of `AHProf` and `RHProf` | m |

Both are constant in time; notebook 02 asserts this with `np.allclose`.
**Unknown whether these are metres above sea level or above antenna.** This
is the prime suspect for the `T_cloud = +4.7 C at 2937 m` anomaly in the
handoff (Sect. 10). Check `TAlts[0]` against the antenna height before
anything else in the scanning work, because `delta_kappa` depends on it.

### Per-chirp scalars, dims `(Chirp,)`

| Name | Meaning | Units |
|---|---|---|
| `SeqIntTime` | sequence integration time per chirp | s |
| `AvgNum` | number of spectral averages | - |
| `DoppLen` | Doppler FFT length | - |
| `MaxVel` | Nyquist velocity | m/s |
| `RangeRes` | range resolution | m |

### File-level scalars

| Name | Meaning |
|---|---|
| `Freq` | radar frequency, GHz |
| `HPBW` | half-power beamwidth, deg |
| `AntDia`, `AntG`, **`AntSep`** | antenna diameter, gain, **separation** |
| `GPSLat`, `GPSLon` | position |
| `DualPol` | 0 single pol, 1 LDR mode, 2 full |
| `ChirpNum` | number of chirps |
| `ProgNo` | chirp program number |
| `ModelNo`, `SWVersion` | instrument model, software version |
| `SampDur`, `CallInt` | sample duration, call interval |
| `FFTInRng`, `FFTWin`, `NoiseFilt`, `SpkFilEna`, `PhaseCor`, `RelPowCor`, `SupPowLev` | processing settings |

`AntSep` is worth reading — the beam-mismatch section currently derives
separation from `GPSLat`/`GPSLon` by haversine, which gave d = 0.

### Loading: the `wanted` set in notebook 02

`load_rpg` builds a keep-list and passes everything else to
`drop_variables`. **Anything not in `wanted` is silently dropped.** Current
contents:

```python
per_chirp = ["ZE", "MeanVel", "SLv", "SLh"]
wanted  = {f"C{c}{v}" for c in (1,2,3) for v in per_chirp}
wanted |= {f"C{c}Range" for c in (1,2,3)}
wanted |= {"Time", "TProf", "TAlts", "AHProf", "HAlts", "SurfPres", "LWP"}
```

For scanning, add:

```python
per_chirp = ["ZE", "MeanVel", "SLv", "SLh", "SLDR", "SpecWidth"]
wanted |= {"Elv", "Azm"}
```

### Filename glob — a blocker, not a detail

```python
files = sorted(glob.glob(str(Path(folder) / "*ZEN.LV1.nc")))
```

`ZEN` is the RPG zenith scan-type token. Scan files carry a different token.
**This glob will return zero files for PPI data and the notebook will fail
at `open_mfdataset`.** Check the actual filenames in the PPI folder first and
make the token a variable:

```python
SCAN_TOKEN = "ZEN"       # or whatever the PPI files actually use
files = sorted(glob.glob(str(Path(folder) / f"*{SCAN_TOKEN}.LV1.nc")))
```

`load_rpg` is defined **three times** with different bodies: nb02 cell 4
(minimal), nb02 cell 12 (with `wanted`/`drop_variables`), nb04 CLEANED cell
50 (minimal again). Only the nb02 cell 12 version matters. Fix all three or
consolidate into `src/`.

---

# PART C — `aligned_LWP2.nc`, dims `(Time, range)`

| Name | Meaning | Units |
|---|---|---|
| `Ze_W`, `Ze_Ka` | reflectivity on the common time grid | dBZ |
| `MeanVel_W`, `MeanVel_Ka` | mean Doppler velocity | m/s |
| `SLv_W`, `SLv_Ka` | sensitivity limit, vertical pol | mm^6/m^3 linear |
| `SLh_W`, `SLh_Ka` | sensitivity limit, horizontal pol | mm^6/m^3 linear |
| `T_cloud` | `TProf` interpolated onto the radar grid | K |
| `AH_cloud` | `AHProf` interpolated onto the radar grid | g/m3 |
| `P_cloud` | barometric pressure profile | hPa |
| `LWP` | MWR LWP, dims `(Time,)` only | g/m2 |

Coordinates carried on `range`: `TAlt`, `HAlt` (leftovers from the
interpolation, harmless).

### Notebook 02 intermediate names

```python
ds_w, ds_ka              # loaded RPG datasets, W and Ka
ze_w, ze_ka              # combined-chirp reflectivity, dBZ
common_times             # np.intersect1d(ze_w.Time, ze_ka.Time)
ze_w_c, ze_ka_c          # reflectivity on common_times   <- "_c" = common
MeanVel_W, MeanVel_Ka    # combined chirps, full time
MeanVel_W_c, MeanVel_Ka_c
SLv_W, SLv_Ka, SLh_W, SLh_Ka
SLv_W_c, SLv_Ka_c, SLh_W_c, SLh_Ka_c
T_alts, T_prof, T_on_range      # -> named "T_cloud"
AH_alts, AH_prof, AH_on_range   # -> named "AH_cloud"
p_surf, surface_pressure, radar_range, pressure_profile, scale_height
P_cloud
lwp
ds_out                   # -> aligned_LWP2.nc
missing_in_w, missing_in_ka, identical, dt_w, dt_ka
lam, max_vel, avg_num, seq_int, T_chirp     # chirp-timing diagnostics
sep, hpbw_ka, hpbw_w, lat_ka, lon_ka, lat_w, lon_w, hs, Dka, Dw, frac
```

### THE SOUNDING INTERPOLATION IS IN NOTEBOOK 02, NOT 03

The handoff puts the `z = r sin(theta)` fix in notebook 03. It belongs in
**notebook 02**. `T_cloud`, `AH_cloud` and `P_cloud` are all built there by
treating `range` as height:

```python
T_on_range  = T_prof.interp(TAlt=ze_w_c["range"], method="linear")   # nb02
AH_on_range = AH_prof.interp(HAlt=ze_w_c["range"], method="linear")  # nb02
pressure_profile[t, r] = p0 * np.exp(-radar_range[r] / scale_height) # nb02
```

Notebook 03 only *consumes* `T_cloud` (for `kappa_field`) and notebook 02_05
consumes `T_cloud`, `P_cloud`, `AH_cloud` (for `alpha_field`). If the height
fix is applied in 03, `alpha_W` and `alpha_Ka` from 02_05 will still be wrong,
because they were computed in the previous notebook. Fix 02, then re-run
02_05 and 03 in order.

The three interpolations also need `z` as a **2-D (ray, gate)** field for an
RHI, which `interp` will not accept as a scalar coordinate. Use
`xr.DataArray(z, dims=("Time","range"))` and `interp` with a broadcast
coordinate, or `np.interp` in a loop over rays.

---

# PART D — `retrieval_input_4km_alpha_kappa_calculated.nc`, dims `(Time, range)`

| Name | Meaning | Units |
|---|---|---|
| `Ze_W`, `Ze_Ka` | reflectivity, sliced 0–4000 m | dBZ |
| `DWR` | `Ze_Ka - Ze_W` | dB |
| `alpha_W`, `alpha_Ka` | gas attenuation, ITU-R P.676-13 | dB/km |
| `kappa_W`, `kappa_Ka` | liquid attenuation, Liebe 1989 dielectric | dB/km/(g/m3) |
| `cloud_mask` | both channels above a constant threshold | bool |
| `common_mask` | significance + 5x5 neighbour, both channels, plus `range >= 500` | bool |
| `combined_mask` | `common_mask & drizzle_mask` | bool |
| `T_cloud` | temperature | K |
| `LWP` | MWR LWP, dims `(Time,)` | g/m2 |
| `MeanVel_W`, `MeanVel_Ka` | mean Doppler velocity | m/s |
| `SLv_W`, `SLv_Ka` | sensitivity limits, linear | mm^6/m^3 |

`SLv_W` / `SLv_Ka` are in the `_alpha_kappa_calculated` file **only**, not in
`retrieval_input_4km.nc`. Notebook 06 Step 1b checks for them and falls back
to an assumed 20 dB SNR if absent.

## The case trap: `alpha_W` vs `alpha_w`

Notebook 03 carries two complete sets of attenuation coefficients under names
that differ only in capitalisation.

| Python name | What it is | Value |
|---|---|---|
| `alpha_W`, `alpha_Ka` | **computed**, ITU-R P.676-13 via `alpha_field` | 0.008–0.655 / 0.005–0.145 dB/km |
| `kappa_W`, `kappa_Ka` | **computed**, Liebe via `kappa_field` | temperature-dependent |
| `alpha_w`, `alpha_ka` | **CONSTANT PLACEHOLDERS** | `np.full_like(dwr, 0.65)` / `0.14` |
| `kappa_w`, `kappa_ka` | **CONSTANT PLACEHOLDERS** | `np.full_like(dwr, 4.29)` / `0.80` |

Both sets are written to netCDF under the **same** variable names
`alpha_W`/`alpha_Ka`/`kappa_W`/`kappa_Ka`, in different files:

- `retrieval_input_4km.nc` → contains the **placeholders**
- `retrieval_input_4km_alpha_kappa_calculated.nc` → contains the **computed** fields

Reading the wrong file gives a plausible LWC field with a constant
`delta_kappa` of 3.49 and `delta_alpha` of 0.51, and nothing warns you.
The `_alpha_kappa_calculated` file is the correct one. Consider deleting the
placeholder cells outright when you copy the folder for the scanning work.

### Notebook 03 intermediate names

```python
ds, ds_alpha, r
dwr, dwr_2km, dwr_2km_diff, dwr_2km_masked
cloud_dwr, fraction
SL_W, SL_Ka, SLv_W_dBZ, SLv_Ka_dBZ
cloud_mask, cloud_mask_2km
sig_W, sig_Ka                           # Ze notnull & > -41 dBZ
neighbour_mask_W, neighbour_mask_Ka     # 5x5 box, >= 13 neighbours
Ze_W_filt, Ze_Ka_filt
common_mask, common_mask_2km
dv, dv_threshold, dv_2km, dv_vals       # MeanVel_W - MeanVel_Ka
gate_flag, profile_flag, drizzle_mask
min_drizzle_height, max_drizzle_height
combined_mask, combined_mask_2km
ze_ka_avg, ze_w_avg                     # 30 s linear-Ze average
dwr_2km_smoothed                        # DWR from the averaged Ze
gate_threshold, dwr_2km_filt            # after despike -> FEEDS THE RETRIEVALS
n_replaced
alpha_w_2km, alpha_ka_2km, kappa_w_2km, kappa_ka_2km   # placeholders, sliced
alpha_W_2km, alpha_Ka_2km, kappa_W_2km, kappa_Ka_2km   # computed, sliced
meanvel_W, meanvel_Ka
min_neighbours = 13, box = (5, 5)
med, mad, sigma_robust                  # MAD-based dv threshold check
exceed_2km, runs, n_exceed, fraction_exceed, max_dv
thresholds, margins, results_const, results_sl
kept_fraction, sigma_dwr, lwp_noise
ds_retrieval, ds_retrieval_2km, ds_retrieval_2km_alpha_kappa
```

---

# PART E — RETRIEVAL OUTPUTS

### `retrieval_output_zhu.nc` (nb04 CLEANED, cell 45)

| netCDF name | Python name | Meaning |
|---|---|---|
| `lwc_guarded` | `lwc_guarded` | LWC, Zhu method, g/m3 |
| `derivative_guarded` | `derivative_guarded` | dDWR/dh, **dB/km** despite the `dB m-1` attribute |
| `m_eff` | `m_eff` | effective sample size per gate, Zhu Eq. 8 |
| `n_windows` | `n_windows` | number of sliding windows per gate (M) |
| `n_used` | `n_used` | window size actually used per profile, gates |
| `dwr_filtered` | `dwr_2km_filt` | the DWR the retrieval actually saw |

### `retrieval_output_zhu_pad{0,1}.nc` (nb04 CLEANED, cell 55)

| netCDF name | Meaning |
|---|---|
| `lwc_guarded`, `dDWR_dh`, `m_eff`, `n_windows` | as above |
| **`N_used`** | **same quantity as `n_used`, capital N** |
| `dwr_filtered`, `common_mask` | inputs |

**`n_used` and `N_used` are the same quantity written under two different
names in two different files.** Notebook 06 reads `n_used` at Step 1 and
`N_used` at cells 43–44. Pick one spelling before adding notebooks 08/09.

### `retrieval_output_dq.nc` (nb05)

| netCDF name | Python name | Meaning |
|---|---|---|
| `lwc_dq` | `lwc_2` | LWC, integral method, g/m3 |
| `dq` | `dq_2` | cumulative liquid water, g/m3 km; `LWP = DQ_top * 1000` |

### `lwc_uncertainty.nc` (nb06, cell 32)

| netCDF name | Meaning |
|---|---|
| `LWC_uncertainty` | 1-sigma per-gate LWC uncertainty, g/m3 |
| `M_effective` | `m_eff` after the Eq. 7 treatment |
| `delta_kappa` | `kappa_W - kappa_Ka` field |
| `n_used` | window size per profile |
| `window_span_m` | `n_used * gate_spacing`, the `dh` of Zhu Eq. 5 |

### `lwc_uncertainty_pad{0,1}.nc` (nb06, cell 43)

`LWC_uncertainty`, `N_used`.

---

# PART F — `src/radar_lwc/` FUNCTION SIGNATURES

```python
# preprocessing
from radar_lwc.preprocessing.chirps          import combine_chirps
from radar_lwc.preprocessing.get_reflectivity import get_reflectivity
from radar_lwc.preprocessing.clean           import drop_duplicate_times
from radar_lwc.preprocessing.align           import align_to_reference_time

combine_chirps(ds, var_prefix="ZE")         # "ZE" | "MeanVel" | "SLv" | "SLh" | "SLDR"
get_reflectivity(ds)                        # -> dBZ on a merged range axis
drop_duplicate_times(ds)
align_to_reference_time(da_ref=, da_target=, method="nearest")

# retrieval
from radar_lwc.retrieval.gas_attenuation import specific_gas_attenuation, alpha_field
from radar_lwc.retrieval.kappa_liquid    import kappa_liquid, kappa_field
from radar_lwc.retrieval.dielectric      import epsilon_water_liebe1989
from radar_lwc.retrieval.average_ze_time import avg_ze_time
from radar_lwc.retrieval.ddwr_dh         import compute_ddwr_dh
from radar_lwc.retrieval.zhu             import retrieve_lwc_zhu
from radar_lwc.retrieval.dq_retrieval    import retrieve_lwc_dq

alpha_field(P_field=, T_field=, AH_field=, freq_ghz=)   # -> dB/km, ITU-R P.676
kappa_field(T_field=, freq_ghz=)                        # -> dB/km/(g/m3)
kappa_liquid(freq_ghz, T_kelvin)                        # scalar, positional
avg_ze_time(da, window_s=30)                            # linear-space Ze average
compute_ddwr_dh(dwr, cloud_mask, N_default=9)
```

### `retrieve_lwc_zhu` — the notebooks use a FIVE-value return

```python
lwc, ddwr_dh, m_eff, n_windows, n_used = retrieve_lwc_zhu(
    dwr        = ...,   # (Time, range) dB
    kappa_w    = ...,   # (Time, range) dB/km/(g/m3)
    kappa_ka   = ...,
    alpha_w    = ...,   # (Time, range) dB/km
    alpha_ka   = ...,
    cloud_mask = ...,   # (Time, range) bool
    N_default       = 9,
    min_cloud_gates = 3,
    max_ddwr_dh     = 50.0,   # dB/km, guard
    max_window_std  = 3.0,    # dB, guard
    edge_pad        = 1,
)
```

**The copy of `zhu.py` in the project folder returns FOUR values** (no
`n_used`) and the handoff says "returns 4 values now". The notebooks unpack
five and work, so the file on disk in the repo is newer than the copy in the
project folder. Use the repo version as authoritative and refresh the project
copy.

Working guard values in nb04 CLEANED, which differ from the defaults:

```python
N_default=9, min_cloud_gates=3, max_ddwr_dh=30.0, max_window_std=0.6, edge_pad=0
```

`retrieve_lwc_zhu_guarded` is also defined **inline** in nb04 CLEANED cell 17
with a two-value return. It is dead code; the cell below it carries the
comment "DONT USE CODE ABOVE HERE RIGHT NOW ANYMORE". Delete it when copying.

### `retrieve_lwc_dq`

```python
lwc, dq = retrieve_lwc_dq(
    dwr=, kappa_w=, kappa_ka=, alpha_w=, alpha_ka=,
    common_mask=,               # note: NOT "cloud_mask" like the Zhu function
    min_cloud_gates=6,
    fit_degree=2,               # 2 or 3
    kappa_mode="tmean",         # "tmean" or "resolved"
)
```

Uses `scipy.ndimage.label` to find contiguous cloud segments, so it handles
multi-layer profiles. `retrieve_lwc_zhu` does **not** — it takes
`cloud_idx.min()` to `cloud_idx.max()` as one segment. On a scanning ray that
crosses several cloud elements at different ranges this merges them and fits a
quadratic across the gaps. More likely to bite for a slant ray than for a
vertical profile.

### `dq_to_lwp`, defined inline in nb05 cell 33

```python
def dq_to_lwp(dq_da, mask):
    dqm  = dq_da.where(mask)
    top  = dqm.ffill("range").isel(range=-1)
    base = dqm.bfill("range").isel(range=0)
    return ((top - base) * 1000.0).where(mask.any("range"))   # g/m2
```

Endpoint difference, so a constant calibration offset between the radars
cancels. **Units of `dq` are `g/m3 * km` along the `range` axis** — so for a
slant ray this returns the along-beam integral, not LWP. Multiply by
`sin(theta)`.

---

# PART G — NOTEBOOK 06 NAMES

### Assumptions, all in Step 1, all replaceable from the LV1 file

```python
SPEED_OF_LIGHT_M_PER_S = 2.99792458e8
frequency_W_hz  = 94.0e9      ; frequency_Ka_hz  = 35.0e9
wavelength_W_m  = 3.189e-3    ; wavelength_Ka_m  = 8.565e-3
prf_W_hz = prf_Ka_hz = 6250.0                     # ASSUMPTION 1, from Hogan 2005
spectral_width_W_m_per_s = spectral_width_Ka_m_per_s = 0.3   # ASSUMPTION 2 -> C{n}SpecWidth
snr_W_linear = snr_Ka_linear = 10**(20.0/10.0)    # ASSUMPTION 3 -> Ze - 10log10(SLv)
dwell_time_s = 30.0                               # ASSUMPTION 4 -> SeqIntTime when scanning
```

`dwell_time_s = 30.0` is the one line the handoff says to change. Correct.
It feeds only `reflectivity_error_db`, i.e. the theoretical sigma.

### Functions

```python
reflectivity_error_db(wavelength_m, prf_hz, spectral_width_m_per_s,
                      snr_linear, dwell_time_s, n_gates)
lwc_uncertainty_single_window(sigma_dwr_db, delta_kappa, window_span_km)
significance_report(lwc_field, uncertainty_field, mask, label)
sigma_dwr_split(dwr, m, rows)        # -> {"edge": (s, n), "interior": (s, n)}
sigma_lwc_ols(sigma_dwr, N, dr_km, dkappa)
sigma_dwr(rows)                      # redefined several times, watch the scope
quad_penalty(N, thickness)           # see warning below
run(ddwr, wstd) ; diagnose(sweep, values, name, unit) ; panel(axrow, rows, label)
```

### Names

```python
mask_retrieval, dwr_filtered, lwc_zhu, n_windows_zhu, lwc_dq
TR, target, to_grid                       # the reindexing idiom, see Part I
r, sp, gate_spacing_m, gate_spacing_km, uniform
n_used, window_span_km, valid
dwr_masked_values, dwr_gate_differences, sigma_dwr_measured_db
assumed_true_gradient_db_per_km, signal_per_gate_db
reflectivity_error_W_db, reflectivity_error_Ka_db
sigma_dwr_theoretical_db, ratio_measured_to_theoretical
delta_kappa_field
lwc_uncertainty_single_measured, lwc_uncertainty_single_theoretical
lwc_uncertainty_wrong_dh                  # the dh sanity check
m_effective, lwc_uncertainty_final
results_zhu, results_dq
DKAPPA, EDGE_PAD, N_DEF, SPLIT, DR_KM, SIG_DWR, DALPHA
_PEN_N, _PEN_T                            # the quadratic-penalty lookup tables
SCALES, sweep, sweep_ddwr, sweep_std, DDWR_VALUES, STD_VALUES
rows_g, rows_s, base_n
ds_uncertainty, output_path
```

### `quad_penalty` is still in the notebook

The handoff (Sect. 10) says the AI-suggested Monte-Carlo quadratic-penalty
correction "was deleted rather than committed" and must not be reintroduced
without independent checking. It is **present** in the uploaded
`06_error_estimation.ipynb`, cells 42–44, as `quad_penalty`, `_PEN_N`,
`_PEN_T`, and it is what writes `lwc_uncertainty_pad{0,1}.nc`. Either the
handoff is stale or the notebook is a working copy that was never cleaned.
Decide which before the scanning branch inherits it.

---

# PART H — CONSTANTS IN ONE PLACE

```python
# geometry / instrument
HPBW_KA_DEG = 0.84          # from the file attribute HPBW
HPBW_W_DEG  = 0.55
EARTH_RADIUS_M = 6_371_000.0
a_ef = (4.0/3.0) * EARTH_RADIUS_M       # Lamer Eq. (1)

# gate spacings, chirp-dependent, VERIFY PER DATASET
#   zenith 20230930 : 22.4 / 27.0 / 37.7 m
#   "Scan 0.5s" prog: 26.5 / 27.5 m
gate_spacing_m = 22.4                   # mode of diff(range), nb06 computes it

# physics
delta_kappa ~ 3.4   dB/km/(g/m3)
delta_alpha ~ 0.51  dB/km
scale_height = 8000.0 m                 # barometric P_cloud in nb02

# placeholders (do NOT ship these)
alpha_w = 0.65 ; alpha_ka = 0.14 ; kappa_w = 4.29 ; kappa_ka = 0.80

# masking
Ze threshold      = -41 dBZ             # constant, cloud_mask and sig_W/sig_Ka
min_neighbours    = 13 of box (5, 5)
dv_threshold      = 0.15 m/s            # 0.1 in some cells, 0.15 in the live one
min_drizzle_height = 550 m ; max_drizzle_height = 4000 m
mask height floor = 500 m               # combined_mask / common_mask
range slice       = slice(0, 4000)

# retrieval
N_default = 9 ; min_cloud_gates = 3 (Zhu) / 6 (DQ)
max_ddwr_dh = 30.0 dB/km ; max_window_std = 0.6 dB ; edge_pad = 0
gate_threshold = 1.0 dB                 # despike
time averaging = 30 s                   # avg_ze_time(window_s=30)
fit_degree = 2 ; kappa_mode = "tmean"

# measured noise (from the handoff, chirp 2, SeqIntTime = 0.137 s)
sigma_DWR single profile        = 0.805 dB
sigma_DWR after 30 s            = 0.29 dB
SIG_DWR interior (used in nb06) = 0.223 dB
sigma_DWR cloud base/top        = 0.426 dB
LDR_THRESHOLD_DB = -20.0                # nb08/09 drafts, insect filter
SNR_MARGIN_DB    = 3.0                  # nb08/09 drafts
```

---

# PART I — IDIOMS THAT REPEAT

### Reindexing onto the DWR grid — nb04 CLEANED, nb06, nb08/09 drafts

```python
TR      = lambda a: a.transpose("Time", "range")
target  = TR(dwr_2km_filt)
to_grid = lambda a: TR(a).reindex_like(target)
mask45  = to_grid(ds_retrieval["common_mask"]).fillna(False).astype(bool)
```

`reindex_like` matches on **labels**, not position. It silently produces NaN
where the range values do not match exactly. After a `range=slice(0, 4000)`
this is fine; after any change to the range axis it is a trap. `mask45` is
named for nothing in particular — it is just the common mask on the target grid.

### sigma from adjacent-gate differences — appears in 03, 04, 05, 06

```python
D     = dwr.where(mask).values
diffs = D[:, 1:] - D[:, :-1]
sigma = np.nanstd(diffs[np.isfinite(diffs)]) / np.sqrt(2)
```

`/sqrt(2)` because `Var(n_{i+1} - n_i) = 2 sigma^2`. Differences along axis 1
= along `range` = along the beam, so this is geometry-safe.

### Profile-level LWP

```python
dz  = float(dwr["range"].diff("range").mean())
lwp = (lwc.where(mask).fillna(0) * dz).sum(dim="range")
```

Geometry-**unsafe**, see Part J item 3.

---

# PART J — SCANNING AUDIT: WHAT THE HANDOFF MISSES

Read the notebooks line by line. Eight items beyond what the handoff lists.
Ordered by how badly they bite.

**1. The sounding interpolation lives in notebook 02, not 03.** See Part C.
The handoff's `z = r sin(theta)` block is correct physics pointed at the wrong
file. `T_cloud`, `AH_cloud`, `P_cloud` are all built in 02; `alpha_W`/
`alpha_Ka` are built in 02_05 from those; 03 only slices them. Fix 02, then
re-run 02_05 and 03.

**2. The `*ZEN.LV1.nc` glob returns nothing for scan files.** Notebook 01/02
loading is not "low risk" — it is a hard stop. Check the real filenames in the
PPI folder before anything else.

**3. LWP integration uses `dz = mean(diff(range))`.** Notebook 04 does this in
at least five places, notebook 05's `dq_to_lwp` does it implicitly. On a slant
ray that integrates along the beam, so LWP comes out too large by
`1/sin(theta)` — 1.41x at 45 deg. The handoff says notebook 04 needs "no
change". The retrieval core needs none; the LWP block does:

```python
dz_vertical = dr_m * np.sin(np.radians(ELEV_DEG))
```

Separately: the MWR `LWP` is a zenith-pointing measurement. Comparing it to a
slant-path integral is not a like-for-like closure test at all. Either keep
the LWP closure test zenith-only, or interleave zenith segments and compare
those.

**4. `range=slice(0, 4000)` is a height slice applied to slant range.** At 45
deg it covers only 0–2828 m of height. To reach the same cloud, the scan needs
`slice(0, 4000/sin(theta))` = `slice(0, 5657)`. Appears ~20 times in notebook
03 and in 04/05/06.

**5. The 500 m mask floor and the 550/4000 m drizzle window are height
thresholds applied to `range`.**

```python
combined_mask = combined_mask.where(combined_mask["range"] >= 500, False)   # nb03
dv = dv.where((dv["range"] >= 550) & (dv["range"] <= 4000))                # nb03
retrieval_mask = common_mask & (ds_retrieval["range"] >= 1500)             # nb04
```

All need `z >= ...` with `z = r sin(theta)`, or they cut the wrong altitudes.
The 500 m floor at 45 deg cuts at 354 m of height, so more near-surface
clutter survives.

**6. The 5x5 `uniform_filter` neighbour test is the one Time-axis operation
that survives scanning — with two conditions.** In `(azimuth, range)` it is
exactly Lamer's range-azimuth spatial filter, and 5 rays at 0.55 deg is 2.75
deg, physically sensible. But (a) it must run **per sweep**, never on
concatenated scans, or the box straddles a scan boundary and correlates rays
minutes apart; (b) `mode="reflect"` mirrors data at the sector edges, which
invents neighbours for the first and last rays. Use `mode="constant",
cval=0` for scans so edge rays are treated as having no neighbours.

**7. The drizzle filter rejects a whole profile.** `profile_flag =
gate_flag.any(dim="range")` then `drizzle_mask = ~profile_flag`. For a scan
that rejects a whole **ray**. Defensible, but one drizzly gate anywhere now
removes an entire azimuth from the ensemble, and the ensemble is the whole
deliverable. Consider a gate-level or segment-level drizzle mask for scanning,
and report how many rays each scan loses.

**8. The constant `-41 dBZ` threshold interacts badly with slant geometry.**
`cloud_mask`, `sig_W` and `sig_Ka` all use a fixed dBZ threshold. Sensitivity
falls as `r^2`, and at 45 deg the same height sits at 1.41x the range, so a
fixed threshold throws away roughly 3 dB more cloud at the top of the layer
than at zenith. Notebook 03 already has the machinery for the alternative
(`SLv_W_dBZ + margin`, swept over 0–19 dB). Switching the scan branch to the
SL-based mask makes the zenith-vs-scan comparison much cleaner, because `SLv`
is range-dependent and absorbs the `r^2` penalty automatically.

Also worth noting, not a bug: `derivative_guarded` carries
`units = "dB m-1"` while the values are dB/km, and the DataArray is named
`dDWR_dh`. For the scanning branch it is `dDWR/ds`. Rename it there or the
thesis figures will be mislabelled.

---

# PART K — WHAT NOTEBOOKS 08/09 EXPECT AND WHAT IS MISSING

The drafts open `retrieval_input_2km_alpha_kappa_calculated_scan.nc` and read:

| Expected | Status |
|---|---|
| `Ze_W`, `Ze_Ka`, `DWR` | exists |
| `alpha_W`, `alpha_Ka`, `kappa_W`, `kappa_Ka` | exists |
| `T_cloud` | exists |
| `SLv_W`, `SLv_Ka` | exists in the `_alpha_kappa_calculated` file only |
| `common_mask` | exists |
| **`elevation`** | **does not exist** — must be added as `Elv` in nb02 |
| **`azimuth`** | **does not exist** — must be added as `Azm` in nb02 |
| **`LDR_Ka`** | **does not exist** — `combine_chirps(ds_ka, var_prefix="SLDR")`, note the S |

The drafts have a fallback for missing `LDR_Ka` (reject `T > 5 C` and
`DWR < -1 dB`) but no fallback for missing `elevation`/`azimuth` — they will
raise `KeyError`. Add `Elv`/`Azm` in notebook 02, save them into
`aligned_LWP2.nc` as `elevation`/`azimuth`, and carry them through 03.

Also fix in the drafts: the input filename (`_2km_` → `_4km_`) and
`AZIMUTH_IS_METEOROLOGICAL`, which is set `True` on an assumption. Confirm
against the RPG documentation or against a known landmark before trusting the
`(x, y)` plan view.

---

# PART L — REGRESSION TEST

Before touching PPI data, run the refactored notebooks with
`SCAN_TYPE = "zenith"` and confirm nothing moved:

```python
old = xr.open_dataset("../data/processed/20230930/retrieval_output_zhu.nc")   # committed
new = xr.open_dataset("../data/processed/20230930_refactor/retrieval_output_zhu.nc")

for v in ["lwc_guarded", "derivative_guarded", "dwr_filtered", "m_eff",
          "n_windows", "n_used"]:
    a, b = old[v].values, new[v].values
    same_nan = np.array_equal(np.isnan(a), np.isnan(b))
    g = np.isfinite(a) & np.isfinite(b)
    print(f"{v:20s} NaN pattern {'OK ' if same_nan else 'DIFF'} "
          f"max|diff| = {np.abs(a[g]-b[g]).max() if g.any() else 0:.3e}")
```

Check the NaN pattern, not just the values. A mask change moves which gates
are retrieved without changing any retrieved number, and `max|diff|` over the
intersection will happily report zero while the coverage has shifted.

Copy the committed output to a fixed path first so a re-run cannot overwrite
the reference. And remember the `pad0`/`pad1` and `lwc_uncertainty*` files
live outside `{DS}/` — see Part A.
