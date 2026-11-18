"""RF and light: decisive numbers for the array16 'beyond RF' study (lens: RF <-> light).
Run: OMP_NUM_THREADS=2 .env/bin/python experiments/beyond_rf/rf_light.py
Writes rf_light_results.json and rf_light.png next to this file. All inputs are stated here;
values marked EST are engineering estimates, not sourced facts."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
def jv(n, x):  # small-argument Bessel series, adequate for x < 1
    from math import factorial
    return sum((-1) ** k / (factorial(k) * factorial(k + n)) * (x / 2) ** (2 * k + n) for k in range(8))

HERE = os.path.dirname(os.path.abspath(__file__))
c, h, mu0, eps0, kB = 2.998e8, 6.626e-34, 4e-7 * np.pi, 8.854e-12, 1.381e-23
alpha, me_eV, lam_c_bar = 1 / 137.036, 510998.95, 3.8616e-13  # reduced Compton wavelength, m
R = {}

# ---------- 0. the field this platform can make ----------
P_ch = 10 ** (5 / 10) * 1e-3            # +5 dBm per channel
P_tot = 4 * P_ch                        # four AD9361 TX channels (rev B: 2x bladeRF xA4) -> 12.6 mW
G = 10 ** (15.5 / 10)                   # panel realised gain at 5.8 GHz (REDESIGN sec. 4)
def field(P, G, r):
    S = P * G / (4 * np.pi * r ** 2)
    E = np.sqrt(2 * 377 * S)            # peak
    return S, E, E / c
S, E, B = field(P_tot, G, 0.3)
R["array_field_0p3m"] = dict(P_tot_mW=P_tot * 1e3, EIRP_mW=P_tot * G * 1e3, S_W_m2=S, S_mW_cm2=S * 0.1,
                             E_peak_V_m=E, B_peak_T=B)
S1, E1, B1 = field(1.0, G, 0.3)        # with a 1 W amplifier (amateur service only)
R["field_1W_PA_0p3m"] = dict(S_W_m2=S1, S_mW_cm2=S1 * 0.1, E_peak_V_m=E1, B_peak_T=B1)
# distance at which the FCC general-population limit (1 mW/cm^2 = 10 W/m^2 above 1.5 GHz) is met
R["mpe_distance_m"] = {"array_12.6mW_15.5dBi": float(np.sqrt(P_tot * G / (4 * np.pi * 10))),
                       "1W_15.5dBi": float(np.sqrt(1.0 * G / (4 * np.pi * 10)))}

# ---------- 1. vacuum: light-by-light ----------
def sigma_gg_cm2(E_eV):                 # Euler-Heisenberg, unpolarised, CM frame, omega << m
    return 973 / (10125 * np.pi) * alpha ** 4 * (E_eV / me_eV) ** 6 * (lam_c_bar * 100) ** 2
R["qed"] = dict(sigma_532nm_cm2=sigma_gg_cm2(2.331), sigma_1eV_cm2=sigma_gg_cm2(1.0),
                sigma_6GHz_photon_cm2=sigma_gg_cm2(h * 6e9 / 1.602e-19),
                thomson_cm2=6.652e-25,
                dn_vacuum_per_T2=4e-24,
                dn_array_field=4e-24 * B ** 2,
                phase_rad_over_1m_532nm=2 * np.pi * 1.0 * 4e-24 * B ** 2 / 532e-9,
                dn_1W_PA=4e-24 * B1 ** 2,
                schwinger_E_V_m=1.32e18, ratio_E_array_to_schwinger=E / 1.32e18)

# ---------- 2a. acousto-optics ----------
def ao(lam, f, v):                      # separation of 0th and 1st order (external, small angle)
    return lam * f / v
parts = {  # name: (v m/s, centre f, bandwidth, aperture m, Psat W EST)
    "TeO2 slow-shear AOD 80 MHz": (650.0, 80e6, 50e6, 4e-3, 0.5),
    "PbMoO4 AOM (Isomet 1205C) 80 MHz": (3630.0, 80e6, 30e6, 2e-3, 1.0),
    "TeO2 longitudinal AOM 200 MHz": (4200.0, 200e6, 50e6, 0.5e-3, 1.0),
    "TeO2 GHz AOM (Brimrose TEF-1000)": (4200.0, 1000e6, 300e6, 76e-6, 1.0),
    "GaP GHz AOM (Brimrose GPF-1500)": (6320.0, 1500e6, 500e6, 76e-6, 1.0),
}
aot = {}
for name, (v, f0, bw, D, Psat) in parts.items():
    lam = 532e-9 if "GaP" not in name else 650e-9
    tau = D / v
    aot[name] = dict(lambda_nm=lam * 1e9, sep_mrad=ao(lam, f0, v) * 1e3, scan_mrad=ao(lam, bw, v) * 1e3,
                     mrad_per_MHz=ao(lam, 1e6, v) * 1e3, transit_us=tau * 1e6, resolvable_spots=bw * tau,
                     eff_at_plus5dBm=float(np.sin(np.pi / 2 * np.sqrt(P_ch / Psat)) ** 2),
                     eff_at_0p5W=float(np.sin(np.pi / 2 * np.sqrt(min(0.5 / Psat, 1))) ** 2),
                     gain_needed_dB_for_Psat=10 * np.log10(Psat / P_ch))
R["acousto_optic"] = aot
# chirp lens: F = v^2/(lam*alpha)
v, lam = 650.0, 532e-9
R["ao_chirp_lens"] = {f"{bw/1e6:.0f}MHz_in_{T*1e6:.0f}us": v ** 2 / (lam * bw / T)
                      for bw, T in [(50e6, 10e-6), (50e6, 100e-6), (20e6, 1e-3)]}
# heterodyne interferometer read by the SDR: phase noise limited displacement
R["heterodyne"] = dict(nm_per_degree=532 / 2 / 360, shot_noise_rad_rtHz_at_10uW=float(
    np.sqrt(2 * 1.602e-19 / (0.3 * 10e-6))))  # 0.3 A/W Si at 532 nm, sqrt(2e/I)

# ---------- 2b. electro-optic ----------
Vpi = 5.0
Vpk = np.sqrt(2 * 50 * P_ch)
m = np.pi * Vpk / Vpi
R["eom"] = dict(Vpk_at_plus5dBm=Vpk, phase_index_rad=m, first_sideband_dBc=20 * np.log10(jv(1, m) / jv(0, m)),
                bulk_LiNbO3_dn_per_V_m=0.5 * 2.2 ** 3 * 30.8e-12,
                freespace_dn_at_array_field=0.5 * 2.2 ** 3 * 30.8e-12 * E / 28,   # field reduced ~1/eps_r EST
                kerr_nitrobenzene_dn=532e-9 * 4.4e-12 * E ** 2)                   # K = 4.4e-12 m/V^2 textbook EST
# RF over fibre: intrinsic direct-modulation link (EST: 0.1 W/A laser slope, 0.8 A/W detector)
sl, rd, RIN, Idc = 0.1, 0.8, 10 ** (-150 / 10), 2e-3
g = (sl * rd) ** 2
N_rin = RIN * Idc ** 2 * 50 / 4          # W/Hz delivered (matched)
N_shot = 2 * 1.602e-19 * Idc * 50 / 4
kT = kB * 290
NF = 10 * np.log10(1 + (N_rin + N_shot + kT) / (g * kT))
R["rf_over_fibre"] = dict(link_gain_dB=10 * np.log10(g), NF_dB=NF,
                          cascade_NF_after_21dB_LNA_dB=10 * np.log10(10 ** 0.15 + (10 ** (NF / 10) - 1) / 10 ** 2.1),
                          fibre_ps_per_m=1e12 * 1.468 / c, mm_per_ps=1e-12 * c / 1.468 * 1e3,
                          deg_at_6GHz_per_K_10m=40e-12 * 0.01 * 6e9 * 360,   # 40 ps/km/K
                          dispersive_TTD_nm_for_352ps_1km_SMF=352 / 17.0,
                          smallest_TTD_bit_38p5ps_mm=38.5e-12 * c / 1.468 * 1e3)

# ---------- 2c. magneto-optic ----------
V_TGG, L = 134.0, 0.01
R["faraday_TGG_1cm"] = dict(rad_at_array_far_field=V_TGG * B * L, rad_at_10uT_loop=V_TGG * 10e-6 * L,
                            shot_noise_rad_rtHz_1mW=float(1 / (2 * np.sqrt(1e-3 / (h * c / 532e-9)))))

# ---------- 2d. NV centre ----------
gam = 28.03e9                            # Hz/T
def loopB(I, r): return mu0 * I / (2 * r)
I5 = 2 * np.sqrt(2 * P_ch / 50)          # shorted loop on a 50 ohm line: current doubles (EST, electrically small)
Bl = loopB(I5, 1e-3)
rabi = gam * Bl / np.sqrt(2)             # S=1 transition, linear B1 perpendicular to the NV axis
R["nv_drive"] = dict(loop_I_peak_mA=I5 * 1e3, B1_uT_1mm_loop=Bl * 1e6, rabi_kHz=rabi / 1e3,
                     rabi_kHz_array_far_field=gam * B / np.sqrt(2) / 1e3,
                     P_for_1MHz_rabi_dBm=10 * np.log10(P_ch * (1e6 / rabi) ** 2 * 1e3))
# wireless relay: patch (6 dBi) at 2.87 GHz, 0.3 m from the panel (panel gain 8.5 dBi EST at 2.87 GHz)
lam287 = c / 2.87e9
Sw = P_tot * 10 ** 0.85 / (4 * np.pi * 0.3 ** 2)
Pr = Sw * 10 ** 0.6 * lam287 ** 2 / (4 * np.pi)
R["nv_wireless"] = dict(S_W_m2=Sw, P_received_dBm=10 * np.log10(Pr * 1e3),
                        rabi_kHz=gam * loopB(2 * np.sqrt(2 * Pr / 50), 1e-3) / np.sqrt(2) / 1e3)
def eta(C, dnu, P_red):                  # CW ODMR shot-noise limit, Lorentzian
    Rph = P_red / (h * c / 680e-9)
    return 4 / (3 * np.sqrt(3)) * dnu / (gam * C * np.sqrt(Rph))
R["nv_sensitivity_T_rtHz"] = {"hobby C=1% 10MHz 10uW": eta(0.01, 10e6, 10e-6),
                              "good  C=3% 1MHz 100uW": eta(0.03, 1e6, 100e-6),
                              "reported LED teaching build": 1e-6, "Schloss 2018 lab": 50e-12}
R["nv_lines"] = dict(zeeman_MHz_per_mT=28.03, earth_50uT_split_MHz=2 * 28.03 * 0.05,
                     span_at_5mT_MHz=2 * 28.03 * 5, temp_kHz_per_K=-74.2)

# ---------- 2e. Rydberg vs a plain antenna ----------
lam24 = c / 2.4e9
Ae = lam24 ** 2 / (4 * np.pi)
Smin = kT * 10 ** 0.15 / Ae              # NF 1.5 dB, 0 dBi
R["rydberg_vs_antenna"] = dict(antenna_E_rms_nV_cm_rtHz=float(np.sqrt(377 * Smin)) * 1e9 / 100,
                               rydberg_uV_cm_rtHz={"Sedlacek2012": 30, "Kumar2017": 3, "2026 heterodyne": 0.76})

# ---------- 3. light -> RF ----------
R["light_to_rf"] = dict(HeNe_mode_beat_MHz={f"L={L*100:.0f}cm": c / (2 * L) / 1e6 for L in (0.15, 0.23, 0.30, 0.45)},
                        DFB_pm_per_GHz_1550=1e12 * 1550e-9 ** 2 / c * 1e9,
                        DFB_GHz_per_K=12.5, self_het_delay_km_for_1MHz_linewidth=c / 1.468 / (np.pi * 1e6) / 1e3,
                        radiation_pressure_nN_per_W=1e9 / c)

json.dump(R, open(os.path.join(HERE, "rf_light_results.json"), "w"), indent=1, default=float)

# ---------- figure ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
f = np.linspace(40e6, 120e6, 200)
for name, v, lam in [("TeO2 slow shear, 532 nm", 650, 532e-9), ("TeO2 slow shear, 650 nm", 650, 650e-9),
                     ("PbMoO4, 532 nm", 3630, 532e-9), ("TeO2 longitudinal, 532 nm", 4200, 532e-9)]:
    ax[0].plot(f / 1e6, np.degrees(lam * f / v), label=name)
ax[0].set(xlabel="RF drive frequency (MHz)", ylabel="angle between 0th and 1st order (deg)",
          title="Acousto-optic deflection: angle = lambda f / v")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
# ODMR spectrum, 4 NV orientations, bias field 3 mT along a generic direction, 14N hyperfine ignored
Bv = 3e-3 * np.array([0.90, 0.40, 0.15]) / np.linalg.norm([0.90, 0.40, 0.15])
axes = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
fr = np.linspace(2.75e9, 2.99e9, 2000); sp = np.ones_like(fr)
for a in axes:
    for s in (+1, -1):
        f0 = 2.87e9 + s * gam * abs(a @ Bv)       # weak-field approximation
        sp -= 0.01 * (5e6 / 2) ** 2 / ((fr - f0) ** 2 + (5e6 / 2) ** 2)
ax[1].plot(fr / 1e9, sp)
ax[1].axvspan(2.842, 2.898, color="orange", alpha=.2, label="one 56 MHz AD9361 capture")
ax[1].set(xlabel="microwave frequency (GHz)", ylabel="red fluorescence (normalised)",
          title="NV ODMR, 3 mT bias: 8 lines = full field vector")
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
P = np.logspace(-3, 0, 100)
for C, dnu, Pr_, lab in [(0.01, 10e6, 10e-6, "C 1 %, 10 MHz, 10 uW red"), (0.03, 1e6, 100e-6, "C 3 %, 1 MHz, 100 uW red")]:
    ax[2].axhline(eta(C, dnu, Pr_) * 1e9, ls="--", color="C0" if C < 0.02 else "C1", lw=3 if C > 0.02 else 1.5, label="shot-noise limit: " + lab)
ax[2].axhline(1000, color="r", label="LED teaching build, reported ~1 uT/rtHz")
ax[2].axhline(0.05, color="g", label="Schloss 2018 lab vector, 50 pT/rtHz")
ax[2].set(yscale="log", ylabel="sensitivity (nT / sqrt(Hz))", xticks=[], title="NV magnetometer sensitivity ladder")
ax[2].legend(fontsize=7); ax[2].grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "rf_light.png"), dpi=130)
print(json.dumps(R, indent=1, default=float))
