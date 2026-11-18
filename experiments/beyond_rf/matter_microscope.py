"""RF, matter and the microwave microscope: decisive numbers for array16 rev B.

Run: OMP_NUM_THREADS=2 ../../.env/bin/python matter_microscope.py
Writes matter_results.json and four PNG figures next to this file. numpy + matplotlib only.
All models are first-order estimates; each block states its assumptions.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
c, eps0, mu0 = 299792458.0, 8.8541878128e-12, 4e-7 * np.pi
h, muB, kB, qe, me = 6.62607015e-34, 9.2740100783e-24, 1.380649e-23, 1.602176634e-19, 9.1093837e-31
eta0 = np.sqrt(mu0 / eps0)
R = {}

# ---------------------------------------------------------------- 1. forces and heating
def forces(P_tot, f, label):
    lam = c / f
    spot = (lam / 2) ** 2                       # best-case focal spot area (diffraction limit)
    S = P_tot / spot                            # W/m^2, upper bound: all power in one lambda/2 spot
    E_pk = np.sqrt(2 * eta0 * S)
    E_rms2 = E_pk ** 2 / 2
    gradE2 = E_rms2 * 2 * (2 * np.pi / lam)     # standing-wave scale: d(E^2)/dx ~ 2k E^2
    out = dict(P_W=P_tot, f_GHz=f / 1e9, spot_cm2=spot * 1e4, S_mW_cm2=S / 10, E_peak_V_m=E_pk,
               F_radiation_absorb_N=P_tot / c, F_radiation_mirror_N=2 * P_tot / c)
    # dielectrophoresis F = 2 pi r^3 eps_m Re[K] grad(E_rms^2)
    out["F_DEP_1mm_bead_in_air_N"] = 2 * np.pi * (0.5e-3) ** 3 * eps0 * 1.0 * gradE2      # K~1 (water/metal bead)
    out["weight_1mm_water_bead_N"] = 4 / 3 * np.pi * (0.5e-3) ** 3 * 1000 * 9.81
    out["F_DEP_1um_particle_in_water_N"] = 2 * np.pi * (1e-6) ** 3 * 78 * eps0 * 0.48 * gradE2
    out["thermal_force_kT_over_r_1um_N"] = kB * 300 / 1e-6
    out["U_ponderomotive_eV"] = qe ** 2 * E_pk ** 2 / (4 * me * (2 * np.pi * f) ** 2) / qe
    out["heat_1mL_water_K_per_min"] = P_tot / 4.186 * 60        # all power absorbed, no losses
    R.setdefault("forces", {})[label] = out

P16 = 16 * 10 ** (5 / 10) * 1e-3            # sixteen channels at +5 dBm = 50.6 mW
forces(P16, 6e9, "16x+5dBm_6GHz")
forces(P16, 2.45e9, "16x+5dBm_2.45GHz")
forces(1.0, 2.45e9, "1W_amplified_2.45GHz")
# microelectrode DEP: +5 dBm into an open 10 um gap (voltage doubles on an open line)
Vpk = 2 * np.sqrt(2 * 50 * 10 ** (5 / 10) * 1e-3)
gap = 10e-6
Erms = Vpk / np.sqrt(2) / gap
sig_w = 2 * np.pi * 2.45e9 * eps0 * 9.6        # effective conductivity of water at 2.45 GHz
Pvol = sig_w * Erms ** 2 * gap ** 3
R["forces"]["microelectrode_10um_gap_+5dBm"] = dict(
    V_peak=Vpk, E_rms_V_m=Erms,
    F_DEP_1um_radius_bead_in_water_N=2 * np.pi * (1e-6) ** 3 * 78 * eps0 * 0.48 * Erms ** 2 / gap,
    thermal_force_N=kB * 300 / 1e-6, heat_in_gap_volume_W=Pvol,
    steady_temp_rise_K=Pvol / (4 * np.pi * 0.6 * gap))
# safety: far-field distance where the panel's main beam falls to the FCC public limit (1 mW/cm^2 above 1.5 GHz)
for f in (2.45e9, 6e9):
    G = 4 * np.pi * 0.1 * 0.1 / (c / f) ** 2          # 100 mm x 100 mm aperture, 100 % efficiency (upper bound)
    eirp = P16 * G
    R.setdefault("safety", {})[f"{f/1e9:g}GHz"] = dict(gain_dBi=10 * np.log10(G), EIRP_W=eirp,
        EIRP_dBm=10 * np.log10(eirp * 1e3), r_public_MPE_m=np.sqrt(eirp / (4 * np.pi * 10.0)))
R["safety"]["1W_amp_into_6dBi_r_m"] = np.sqrt(4.0 / (4 * np.pi * 10.0))

# ---------------------------------------------------------------- 2. spin resonance
g_e, g_dpph = 2.002319, 2.0036
fs = np.array([0.4, 0.915, 1.42, 2.45, 2.87, 3.0, 5.8, 6.0]) * 1e9
R["esr"] = dict(field_mT={f"{f/1e9:g}": h * f / (g_dpph * muB) * 1e3 for f in fs},
                MHz_per_mT=g_dpph * muB / h * 1e-9,
                proton_field_T_for={f"{f/1e9:g}GHz": f / 42.577478e6 for f in (0.4, 0.047e9 / 1e9, 6.0)
                                    for f in [f * 1e9]})
# DPPH signal: Curie susceptibility, Lorentzian line, resonator perturbation
N = 6.02214076e23 / 394.32 * 1.4e6               # spins per m^3 (density about 1.4 g/cm^3: estimate)
chi0 = mu0 * N * g_dpph ** 2 * muB ** 2 * 0.75 / (3 * kB * 295)
dB_pp = 0.2e-3                                    # 2 G peak-to-peak (literature range 1.5-4.7 G)
dB_hwhm = dB_pp * np.sqrt(3) / 2
sig = {}
for f in (0.4e9, 0.915e9, 2.45e9, 5.8e9):
    B0 = h * f / (g_dpph * muB)
    chi2 = chi0 * B0 / (2 * dB_hwhm)              # chi'' at line centre
    row = dict(B0_mT=B0 * 1e3, chi_pp=chi2, linewidth_MHz_pp=dB_pp * g_dpph * muB / h / 1e6)
    for name, eta_f, Q in (("loop_gap_Q200_fill0.05", 0.05, 200), ("microstrip_Q50_fill0.005", 0.005, 50),
                           ("shorted_line_noQ_fill0.005", 0.005, 1)):
        dQinv = eta_f * chi2
        # critically coupled reflection resonator: Gamma changes by about (Q/2) * delta(1/Q) for small changes
        dGamma = min(1.0, 0.5 * Q * dQinv)
        row[name] = dict(delta_inv_Q=dQinv, delta_reflection=dGamma, dB_down=20 * np.log10(dGamma))
    sig[f"{f/1e9:g}GHz"] = row
R["esr"]["dpph"] = dict(chi0=chi0, spins_per_mg=6.02214076e23 / 394.32 * 1e-3, signal=sig)
Ms = 0.175                                        # mu0*Ms of YIG in tesla (about 1750 G)
R["yig"] = dict(GHz_per_T=28.025, f_min_sphere_GHz=28.025 * Ms / 3,
                field_mT_for={f"{f:g}GHz": f / 28.025 * 1e3 for f in (1.7, 2.45, 4.0, 6.0)},
                Q_at_3GHz_for_dH_0p5_Oe=3e9 / (2.8025e6 * 0.5), Q_at_3GHz_for_dH_1p5_Oe=3e9 / (2.8025e6 * 1.5))
R["nv"] = dict(zero_field_GHz=2.870, MHz_per_mT=28.03, field_mT_to_reach_6GHz=(6.0 - 2.87) / 0.02803)

# ---------------------------------------------------------------- 3a. dielectric spectroscopy
def water(f, T=25.0):
    es, einf, tau = 78.4, 5.2, 8.27e-12          # single Debye, 25 C
    return einf + (es - einf) / (1 + 1j * 2 * np.pi * f * tau)   # returns eps' - j eps''
def pen_depth(er, f):                            # er = eps' - j eps'' ; 1/e field depth
    k = 2 * np.pi * f / c * np.sqrt(er)
    return 1 / abs(k.imag)
mats = {"PTFE": 2.1 * (1 - 0.0003j), "dry concrete": 4.5 * (1 - 0.011j), "ice": 3.2 * (1 - 5e-4j),
        "dry wood (mid)": 2.0 * (1 - 0.03j), "FR-4": 4.4 * (1 - 0.02j)}
ft = np.array([0.4, 0.915, 2.45, 5.8]) * 1e9
tab = {}
for f in ft:
    ew = water(f)
    row = {"water": dict(eps_r=ew.real, eps_i=-ew.imag, tan_d=-ew.imag / ew.real,
                          depth_mm=pen_depth(ew, f) * 1e3,
                          reflect_from_air=abs((1 - np.sqrt(ew)) / (1 + np.sqrt(ew))))}
    for m, e in mats.items():
        row[m] = dict(eps_r=e.real, depth_m=pen_depth(e, f), reflect_from_air=abs((1 - np.sqrt(e)) / (1 + np.sqrt(e))))
    # moist soil/wood by the complex refractive index mixing rule (CRIM): sqrt(e) = sum of volume fractions
    for th in (0.05, 0.15, 0.30):
        por = 0.45
        e = ((1 - por) * np.sqrt(4.0) + th * np.sqrt(ew) + (por - th) * 1.0) ** 2
        topp = np.roots([4.3e-6, -5.5e-4, 2.92e-2, -5.3e-2 - th])
        topp = float(topp[np.isreal(topp)].real.min())
        row[f"soil theta={th}"] = dict(eps_r=e.real, eps_i=-e.imag, topp_eps=topp,
                                       loss_dB_per_m_one_way=8.686 / pen_depth(e, f), depth_m=pen_depth(e, f))
    tab[f"{f/1e9:g}GHz"] = row
R["dielectric"] = tab

# ---------------------------------------------------------------- 3b. resonant probes
# (i) Lumped split-ring probe: f = f0 / sqrt(1 + q (eps - 1)); q = share of the capacitance that the sample fills.
# Calibrated on Isakov et al. 2017 (Sci. Rep.): 1.933 GHz in air -> 1.731 GHz on eps = 2.5.
q = ((1.933 / 1.731) ** 2 - 1) / 1.5
eps = np.linspace(1, 30, 400)
f_srr = 1.933 / np.sqrt(1 + q * (eps - 1))
def srr_sens(e):                                 # MHz per unit eps
    return 1933 * q / 2 / (1 + q * (e - 1)) ** 1.5
def freq_noise(f_MHz, Q, snr_dB):                # resonance-fit precision ~ f / (2 Q SNR_amplitude)
    return f_MHz / (2 * Q * 10 ** (snr_dB / 20))
R["srr_probe"] = dict(q=q, sens_MHz_per_eps={f"eps={e}": srr_sens(e) for e in (1, 2.5, 5, 10, 25)},
    freq_precision_kHz_Q50_SNR40dB=freq_noise(1933, 50, 40) * 1e3,
    d_eps_noise_limit_at_eps3=freq_noise(1933, 50, 40) / srr_sens(3),
    d_eps_from_0p1pct_freq_repeatability_at_eps3=1.9 / srr_sens(3),
    d_eps_at_eps25_same=1.9 / srr_sens(25))
# small-sample cavity perturbation (rod along E, no depolarisation): df/f = -(eps'-1) Vs/(2 Veff), d(1/Q) = eps'' Vs/Veff
Veff = 20e-3 * 20e-3 * 1.6e-3 / 2                # 20 mm microstrip-ring-class effective mode volume (estimate)
R["perturbation"] = {f"sample {d} mm cube": dict(Vs_over_Veff=(d * 1e-3) ** 3 / Veff,
        df_MHz_at_2p45_eps3=-2450 * (3 - 1) * (d * 1e-3) ** 3 / (2 * Veff),
        df_MHz_at_2p45_water=-2450 * 3 * (77 - 1) / (77 + 2) * (d * 1e-3) ** 3 / (2 * Veff),  # sphere-depolarised
        dQinv_water=9.6 * 9 / (77 + 2) ** 2 * (d * 1e-3) ** 3 / Veff) for d in (0.1, 0.3, 1.0)}

# (ii) Sharp tip on a quarter-wave coax resonator: image-charge series (Gao and Xiang, Rev. Sci. Instrum. 1998)
def tip_series(er, g_over_R, nmax=4000):
    b = (er - 1) / (er + 1)
    a1 = 1 + g_over_R
    a, t, s = a1, 1.0, 0.0
    for _ in range(nmax):
        s += b * t / (a1 + a)
        t = b * t / (1 + g_over_R + a)
        a = 1 + g_over_R - 1 / (1 + g_over_R + a)
        if t < 1e-14: break
    return s                                      # df/f = -A * s
def A_const(R0, f):                               # A = 4 pi eps0 R0 / (2 C_eff), C_eff = half the line capacitance
    l = c / f / 4
    return 4 * np.pi * eps0 * R0 / (2 * 0.5 * l / (50 * c))
tip = {}
for R0 in (1e-6, 10e-6, 50e-6, 150e-6):
    A = A_const(R0, 2.45e9)
    row = dict(A=A)
    for er in (2.1, 4.4, 10, 77):
        s0 = tip_series(er, 0.0)
        row[f"eps={er}"] = dict(df_kHz_contact=-A * s0 * 2.45e6,
            df_kHz_gap_0p1R=-A * tip_series(er, 0.1) * 2.45e6, df_kHz_gap_1R=-A * tip_series(er, 1.0) * 2.45e6)
    # smallest eps step near eps=4.4 for a frequency precision of f/(2 Q SNR): Q = 150, SNR 60 dB
    df_min = 2.45e6 / (2 * 150 * 1000)            # kHz
    de = 0.05
    slope = A * (tip_series(4.4 + de, 0) - tip_series(4.4 - de, 0)) / (2 * de) * 2.45e6
    row["d_eps_noise_limit_near_4p4"] = df_min / slope
    # gap error that mimics a 5 % change in eps at eps=4.4, in contact region
    dg = 0.01
    gslope = A * (tip_series(4.4, 0.0) - tip_series(4.4, dg)) / dg * 2.45e6    # kHz per (gap/R0)
    row["gap_change_um_equal_to_5pct_eps"] = 0.05 * 4.4 * slope / gslope * R0 * 1e6
    tip[f"R0={R0*1e6:g}um"] = row
R["tip_probe"] = tip
R["tip_probe"]["freq_precision_kHz_Q150_SNR60dB"] = 2.45e6 / (2 * 150 * 1000)

# ---------------------------------------------------------------- 3c. array imaging point-spread functions
def psf(tx, rx, freqs, tgt, axis, span, n=401, pairs=None):
    k = 2 * np.pi * freqs / c
    if pairs is None:
        it, ir = np.meshgrid(np.arange(len(tx)), np.arange(len(rx)), indexing="ij")
        it, ir = it.ravel(), ir.ravel()
    else:
        it, ir = pairs
    T, Rx = tx[it], rx[ir]
    d0 = np.linalg.norm(T - tgt, axis=1) + np.linalg.norm(Rx - tgt, axis=1)
    s = np.exp(-1j * np.outer(k, d0))             # (nf, npairs), unit amplitudes
    u = np.linspace(-span, span, n)
    img = np.empty(n)
    for i, ui in enumerate(u):
        p = tgt.copy(); p[axis] += ui
        d = np.linalg.norm(T - p, axis=1) + np.linalg.norm(Rx - p, axis=1)
        img[i] = abs(np.sum(s * np.exp(1j * np.outer(k, d))))
    img /= img.max()
    above = np.where(img >= 0.5)[0]               # -6 dB two-way amplitude = -3 dB per one-way pass
    j = n // 2; lo = j; hi = j
    while lo > 0 and img[lo - 1] >= 0.5: lo -= 1
    while hi < n - 1 and img[hi + 1] >= 0.5: hi += 1
    width = u[hi] - u[lo]
    side = img[np.r_[0:lo, hi + 1:n]].max() if (lo > 0 or hi < n - 1) else 0.0
    return u, img, width, side
g4 = (np.arange(4) - 1.5) * 0.025
panel = np.array([[x, y, 0.0] for x in g4 for y in g4])
quads = np.array([[x, y, 0.0] for x in (-0.025, 0.025) for y in (-0.025, 0.025)])
wide = np.linspace(2.3e9, 6.0e9, 38)
img_res, curves = {}, {}
for z0 in (0.05, 0.15, 0.5):
    tgt = np.array([0.0, 0.0, z0])
    for name, tx, rx, fr in (("16x16 MIMO, 2.3-6 GHz", panel, panel, wide),
                             ("16x16 MIMO, 6 GHz only", panel, panel, np.array([6e9])),
                             ("16x16 MIMO, 2.45 GHz only", panel, panel, np.array([2.45e9])),
                             ("4 quadrants x 4 quadrants, 2.3-6 GHz", quads, quads, wide)):
        u, im, wx, sx = psf(tx, rx, fr, tgt, 0, 0.3)
        row = dict(cross_range_mm=wx * 1e3)
        if len(fr) > 1:
            uz, imz, wz, sz = psf(tx, rx, fr, tgt, 2, 0.04 if z0 < 0.1 else 0.3)
            row["range_mm"] = wz * 1e3
        img_res[f"{name} @ {z0*1e3:g} mm"] = row
        if z0 == 0.15: curves[name] = (u, im)
# synthetic aperture: one transmit/receive element carried over 200 x 200 mm on the printer bed, 10 mm steps
gs = np.arange(-0.1, 0.1001, 0.01)
sar = np.array([[x, y, 0.0] for x in gs for y in gs])
idx = np.arange(len(sar))
for z0 in (0.05, 0.15):
    tgt = np.array([0.0, 0.0, z0])
    for name, fr in (("SAR 200 mm, 2.3-6 GHz", wide), ("SAR 200 mm, 6 GHz only", np.array([6e9]))):
        u, im, wx, sx = psf(sar, sar, fr, tgt, 0, 0.1, pairs=(idx, idx))
        img_res[f"{name} @ {z0*1e3:g} mm"] = dict(cross_range_mm=wx * 1e3)
        if z0 == 0.15: curves[name] = (u, im)
R["imaging"] = img_res
R["imaging_rules"] = dict(range_res_mm_air_3p7GHz=c / (2 * 3.7e9) * 1e3, range_res_mm_56MHz=c / (2 * 56e6) * 1e3,
    half_wavelength_mm={f"{f:g}GHz": c / f / 2e6 for f in (2.3, 2.45, 5.8, 6.0)})
# ground-penetrating depth: how deep a 70 dB two-way soil-loss allowance reaches (the limit in practice is the
# direct antenna-to-antenna leak and surface clutter, not thermal noise). Loss = water relaxation (CRIM mix)
# plus ionic conductivity sigma: alpha = 8.686 * sigma * eta0 / (2 sqrt(eps)) dB/m. Estimates.
bud = {}
for fG in (0.4, 0.915, 2.45, 5.8):
    for th, sg in ((0.05, 0.002), (0.15, 0.02), (0.30, 0.05)):
        f = fG * 1e9
        e = ((1 - 0.45) * 2.0 + th * np.sqrt(water(f)) + (0.45 - th)) ** 2
        aL = 8.686 / pen_depth(e, f) + 8.686 * sg * eta0 / (2 * np.sqrt(e.real))
        bud[f"{fG} GHz, theta={th}, sigma={sg} S/m"] = dict(one_way_dB_per_m=aL, depth_m_for_70dB_two_way=70 / (2 * aL),
            range_res_mm_1GHz_bandwidth=c / (2 * 1e9 * np.sqrt(e.real)) * 1e3)
R["gpr_budget"] = bud

# ---------------------------------------------------------------- figures
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ff = np.logspace(8, 10.5, 300); ew = water(ff)
ax[0].semilogx(ff / 1e9, ew.real, label="eps'"); ax[0].semilogx(ff / 1e9, -ew.imag, label="eps''")
ax[0].axvspan(0.4, 6, alpha=0.15, color="g"); ax[0].set_xlabel("GHz"); ax[0].set_title("Water at 25 C (Debye); band shaded"); ax[0].legend()
ax[1].loglog(ff / 1e9, [pen_depth(e, f) * 1e3 for e, f in zip(ew, ff)]); ax[1].axvspan(0.4, 6, alpha=0.15, color="g")
ax[1].set_xlabel("GHz"); ax[1].set_ylabel("mm"); ax[1].set_title("1/e field penetration depth in water")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "water_debye.png"), dpi=110); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].plot(eps, f_srr); ax[0].set_xlabel("sample permittivity"); ax[0].set_ylabel("GHz")
ax[0].set_title(f"Split-ring probe (q = {q:.3f}, fitted to Isakov 2017)")
ee = np.linspace(1.5, 80, 60)
for gr in (0, 0.1, 0.5, 1, 3):
    ax[1].plot(ee, [-tip_series(e, gr) for e in ee], label=f"gap = {gr} R0")
ax[1].set_xlabel("sample permittivity"); ax[1].set_ylabel("df / (A f)"); ax[1].set_title("Tip probe: shift against permittivity"); ax[1].legend()
gg = np.logspace(-2, 1.3, 60)
for er in (2.1, 4.4, 10, 77):
    ax[2].semilogx(gg, [tip_series(er, g) / tip_series(er, 0) for g in gg], label=f"eps = {er}")
ax[2].set_xlabel("gap / tip radius"); ax[2].set_ylabel("signal / contact signal"); ax[2].set_title("Tip probe: lift-off"); ax[2].legend()
fig.tight_layout(); fig.savefig(os.path.join(HERE, "probe_shift.png"), dpi=110); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
fa = np.linspace(0.4, 6, 100)
ax[0].plot(fa, h * fa * 1e9 / (g_dpph * muB) * 1e3, label="electron spin (g = 2.0036)")
ax[0].plot(fa[fa > 28.025 * Ms / 3], fa[fa > 28.025 * Ms / 3] / 28.025 * 1e3, "--", label="YIG sphere (above saturation)")
ax[0].plot(fa[fa > 2.87], (fa[fa > 2.87] - 2.87) / 0.02803, ":", label="NV centre, upper branch")
ax[0].set_xlabel("GHz"); ax[0].set_ylabel("field, mT"); ax[0].set_title("Resonance field across the band"); ax[0].legend()
fx = np.linspace(-15, 15, 600); wid = dB_hwhm * g_dpph * muB / h / 1e6
ax[1].plot(fx, -2 * fx / wid / (1 + (fx / wid) ** 2) ** 2 / wid)
ax[1].set_xlabel("offset from resonance, MHz"); ax[1].set_title("DPPH derivative line in a frequency sweep (0.2 mT wide)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "spin_resonance.png"), dpi=110); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4))
for nme, (u, im) in curves.items():
    ax.plot(u * 1e3, 20 * np.log10(np.maximum(im, 1e-3)), label=nme)
ax.set_ylim(-30, 1); ax.set_xlabel("cross-range, mm"); ax.set_ylabel("dB"); ax.set_title("Point response, target 150 mm in front")
ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(HERE, "imaging_psf.png"), dpi=110); plt.close(fig)

def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (np.floating, float)): return float(f"{float(o):.4g}")
    if isinstance(o, (np.integer,)): return int(o)
    return o
json.dump(clean(R), open(os.path.join(HERE, "matter_results.json"), "w"), indent=1)
print(json.dumps(clean(R), indent=1))
