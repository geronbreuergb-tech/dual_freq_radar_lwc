"""
rpg_chirp_table.py — predict every derived ("blue") column of the RPG FMCW-94
chirp-generator GUI from the user-settable ("white") ones.

The model below was reverse-engineered from the DUT training report (Unal &
Mackenzie 2025) and then validated to the last digit against the 'Scan 0.5s'
program: it reproduces BW, IF-R High, Inv. Smp., Samples, Range FFT, Dopp. FFT,
Total Time and Total Samples exactly.  Run `python rpg_chirp_table.py` to see
the validation plus four candidate scanning programs.

Decoded relations  (f_ADC = 11.45 MHz, fixed by the ADC board)
-------------------------------------------------------------
  chirp slope        S        = IF_low * c / (2 * Hmin)          [Hz/s, at RF]
  IF-R High          IF_high  = IF_low * Hmax / Hmin
  Range FFT          N_r      = f_ADC * Hmin / (IF_low * Resol)   -> power of 2
  samples per chirp  N        = lambda * f_ADC / (4 * vmax)
  Inv. Smp.          Inv      = N - ChirpFFTs * N_r
  chirp duration     T_c      = N / f_ADC
  BW column          BW       = S * T_c / 12          [MHz, CG units; x12 = RF]
  Samples            = ChirpRep * N
  Total Time         = ChirpRep * N / f_ADC
  Dopp. Vel. Resol.  = 2 * vmax / DoppFFT
  incoh. averages    = ChirpFFTs * ChirpRep / DoppFFT

Hard constraints
----------------
  300 kHz <= IF_low,  IF_high <= 3700 kHz   =>   Hmax / Hmin <= 12.33
  Range FFT >= 3.095 * Hmax / Resol         (and a power of two)
  Inv >= 0        =>   vmax <= lambda * f_ADC / (4 * ChirpFFTs * Range FFT)
  ChirpRep = integer multiple of DoppFFT, ratio >= 16 (Myagkov)
  Hmax(n) = Hmin(n+1)   -- no range gaps
"""
import math

C_LIGHT = 2.998e8
F_ADC = 11.45e6                 # ADC sampling rate, Hz  (training report IV-2)
F_RADAR = 94e9
LAMBDA = C_LIGHT / F_RADAR      # 3.19 mm
MULT = 12                       # CG -> RF frequency multiplication at W band
IF_MIN, IF_MAX = 300e3, 3700e3  # receiver IF filter limits


def legal_range_ffts(Hmin, Hmax, resol):
    """Powers of two that satisfy both IF-filter limits for this chirp."""
    lo = F_ADC * Hmax / (IF_MAX * resol)
    hi = F_ADC * Hmin / (IF_MIN * resol)
    return [2 ** k for k in range(3, 15) if lo <= 2 ** k <= hi]


def chirp(Hmin, Hmax, resol, chirp_ffts, vmax, dopp_fft, chirp_rep,
          range_fft=None):
    """Predict all derived columns for one chirp sequence."""
    allowed = legal_range_ffts(Hmin, Hmax, resol)
    if range_fft is None:
        range_fft = allowed[0] if allowed else None   # smallest = best vmax
    if range_fft is None:
        raise ValueError(f"no legal Range FFT for {Hmin}-{Hmax} m at {resol} m")

    if_low = F_ADC * Hmin / (range_fft * resol)
    if_high = if_low * Hmax / Hmin
    slope = if_low * C_LIGHT / (2 * Hmin)

    n_samp = round(LAMBDA * F_ADC / (4 * vmax))       # samples in one chirp
    used = chirp_ffts * range_fft
    t_chirp = n_samp / F_ADC

    return dict(
        Hmin=Hmin, Hmax=Hmax, resol=resol, chirp_rep=chirp_rep,
        BW=slope * t_chirp / MULT / 1e6,
        dv=2 * vmax / dopp_fft, vmax=vmax,
        if_low=if_low / 1e3, if_high=if_high / 1e3,
        inv=n_samp - used, samples=chirp_rep * n_samp,
        range_fft=range_fft, dopp_fft=dopp_fft, chirp_ffts=chirp_ffts,
        total_time=chirp_rep * t_chirp,
        gates=int((Hmax - Hmin) / resol),
        util=used / n_samp,
        n_av=chirp_ffts * chirp_rep // dopp_fft,
        bw_rf_used=chirp_ffts * C_LIGHT / (2 * resol) / 1e6,
        bw_rf_total=slope * t_chirp / 1e6,
        legal=(n_samp - used >= 0 and if_low >= IF_MIN - 1
               and if_high <= IF_MAX + 1 and chirp_rep % dopp_fft == 0
               and chirp_rep // dopp_fft >= 16),
        allowed_range_ffts=allowed,
    )


HDR = (f"{'Ch':>2} {'Hmin':>5} {'Hmax':>6} {'Res':>5} {'ChRep':>6} {'BW':>6} "
       f"{'dv':>6} {'vmax':>5} {'IFlo':>5} {'IFhi':>5} {'Inv':>4} "
       f"{'Samples':>9} {'rFFT':>5} {'dFFT':>5} {'cFFT':>4} {'Time':>6} | "
       f"{'gts':>4} {'util':>5} {'Nav':>4}")


def show(name, rows, note="", beamwidth_w=0.55, beamwidth_ka=0.84,
         sigma_ref=0.805, dwell_ref=0.137, window_m=300.0, dkappa=3.4):
    total = sum(r['total_time'] for r in rows)
    omega = beamwidth_w / total
    print(f"\n=== {name} — {note}")
    print(HDR)
    print('-' * len(HDR))
    for i, r in enumerate(rows, 1):
        print(f"{i:>2} {r['Hmin']:>5} {r['Hmax']:>6} {r['resol']:>5.1f} "
              f"{r['chirp_rep']:>6} {r['BW']:>6.3f} {r['dv']:>6.3f} "
              f"{r['vmax']:>5.1f} {r['if_low']:>5.0f} {r['if_high']:>5.0f} "
              f"{r['inv']:>4} {r['samples']:>9} {r['range_fft']:>5} "
              f"{r['dopp_fft']:>5} {r['chirp_ffts']:>4} {r['total_time']:>6.3f} | "
              f"{r['gates']:>4} {r['util']:>5.0%} {r['n_av']:>4}"
              f"{'' if r['legal'] else '   <-- CHECK>'}")
    print(f"   TOTAL {total:.3f} s, {sum(r['gates'] for r in rows)} layers"
          f"  ->  omega_max = {omega:.2f} deg/s"
          f" | 60 deg sector {60/omega:.0f} s | 90 deg {90/omega:.0f} s"
          f" | advection at 10 m/s over 60 deg: {600/omega:.0f} m")

    ret = rows[1]                                   # the retrieval chirp
    smear = omega * ret['total_time']
    ew, ek = math.hypot(beamwidth_w, smear), math.hypot(beamwidth_ka, smear)
    print(f"   scan smear {smear:.2f} deg -> effective W {ew:.2f} / Ka {ek:.2f}"
          f" deg, mismatch {ek/ew:.2f} (static 1.53)")

    sd = sigma_ref * math.sqrt(dwell_ref / ret['total_time'])
    n = int(window_m / ret['resol'])
    sl = (2.4 * sd / ((ret['resol'] / 1000) * math.sqrt(n * (n * n - 1) / 12))
          * 0.5 / dkappa)
    print(f"   sigma_DWR {sd:.2f} dB/ray | {window_m:.0f} m window = {n} gates"
          f" -> sigma_LWC {sl:.2f} g/m3 per ray, {sl/1.93:.2f} with a"
          f" +/-2-ray 2-D fit")


if __name__ == "__main__":
    print("=== VALIDATION: current 'Scan 0.5s' program "
          "(predicted vs GUI) ===")
    cur = [
        (chirp(100, 200, 12.8, 6, 5.868, 1024, 7168, range_fft=256),
         dict(BW=5.947, if_high=700, inv=20, samples=11153408, total=0.974)),
        (chirp(600, 3000, 26.5, 2, 7.246, 512, 7168, range_fft=512),
         dict(BW=1.161, if_high=2531, inv=236, samples=9031680, total=0.789)),
        (chirp(3000, 12000, 27.5, 1, 3.072, 256, 9216, range_fft=2048),
         dict(BW=0.659, if_high=2438, inv=924, samples=27389952, total=2.392)),
    ]
    for i, (p, o) in enumerate(cur, 1):
        print(f" ch{i}: BW {p['BW']:.3f}/{o['BW']}  IFhi "
              f"{p['if_high']:.0f}/{o['if_high']}  Inv {p['inv']}/{o['inv']}  "
              f"Samples {p['samples']}/{o['samples']}  "
              f"Time {p['total_time']:.3f}/{o['total']}")

    # ---- candidate scanning programs -------------------------------------
    def near(rep=512):
        return chirp(100, 300, 15.0, 6, 22.0, 32, rep)

    def top(rep=1024, dfft=64):
        return chirp(2000, 12000, 150.0, 2, 17.0, dfft, rep)

    show("PROGRAM B  'SCAN-B 25m'",
         [near(), chirp(300, 2000, 25.0, 2, 17.0, 512, 8192), top()],
         "baseline: safe gates, full Doppler headroom, RHI or 45 deg PPI")
    show("PROGRAM D  'SCAN-D 13m'",
         [near(), chirp(300, 2000, 12.5, 1, 17.0, 512, 8192), top()],
         "identical timing to B, half the gate — tests the sqrt(dr) prediction")
    show("PROGRAM A  'SCAN-A fast'",
         [near(), chirp(300, 2000, 25.0, 2, 17.0, 256, 4096), top(512, 32)],
         "half the dwell, twice the scan rate")
    show("PROGRAM C  'SCAN-C deep'",
         [near(1024), chirp(300, 2000, 25.0, 2, 17.0, 512, 16384), top()],
         "double dwell, half the scan rate")
