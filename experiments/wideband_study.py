"""Wideband study: what limits the array's frequency range, and what a KrakenSDR-like design trades.

    python experiments/wideband_study.py

Writes experiments/wideband_results.json and the figures in experiments/wideband/.
The simulator's RF parts are frequency-flat, so its runs off 5.8 GHz show the array as if the
LNAs, combiners and antennas were wideband: what remains is the architecture (delay-line phase
states, aperture, element pitch) and the signal processing.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, ".")
from array16 import Array, antenna as A, board as B, dsp  # noqa: E402
from array16.control import H4, H16, phase_states, snapshot_words, state_lag, toggle_factor  # noqa: E402
from array16.sim import SimBackend, Source  # noqa: E402

OUT = Path("experiments/wideband")
F0 = B.F0
Z0 = 50.0


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-30))


# ---- 1. the 180-degree line away from 5.8 GHz: what the snapshot decode pays in noise -------------
def decode_noise_gain(f):
    """Noise power after decoding, relative to exact 0/180 switching (mean over elements)."""
    t = toggle_factor(0, f)
    out = {}
    for name, H in (("hackrf_16", H16), ("quad_4", H4)):
        M = np.where(H < 0, t[: H.shape[1]][None, :], 1)
        ref = np.trace(np.linalg.inv(H.T @ H)).real
        out[name] = float(np.trace(np.linalg.inv(M.conj().T @ M)).real / ref)
    return out


# ---- 2. direction finding against frequency, in the simulator -------------------------------------
def snapshot_assuming_180(a):
    """What the code did before it knew about the delay lines: decode as if every toggle were -1."""
    words, _ = snapshot_words(0, a.hackrf, False)
    seg, idx = a._slots(words)
    tone = dsp.estimate_tone(seg, idx, a.b.fs)
    return dsp.decode_coherent(dsp.slot_means(seg, idx, a.b.fs, tone).mean(0), a.hackrf, F0, 0)


def df_vs_frequency(freqs, kind, directions, seeds=(21, 22), compare_old=()):
    rows = []
    for f in freqs:
        new, old = [], []
        for seed in seeds:
            b = SimBackend(kind, seed=seed)
            a = Array(b)
            a.tune(f)
            b.sources = [Source(0, 0, -75)]
            a.calibrate()
            ref_old = snapshot_assuming_180(a) if f in compare_old else None
            for d in directions:
                b.sources = [Source(*d, power_dbm=-75)]
                new.append(dsp.angle_between(a.locate(), d))
                if ref_old is not None:
                    z = snapshot_assuming_180(a) / ref_old
                    z = z if a.hackrf else z / a.cal.chan[B.QUAD]
                    old.append(dsp.angle_between(B.angles(*dsp.locate(z, f)), d))
        st = lambda e: dict(median=round(float(np.median(e)), 2), p95=round(float(np.percentile(e, 95)), 2))
        row = dict(f_GHz=f / 1e9, frequency_aware=st(new))
        if old:
            row["decode_assuming_180"] = st(old)
        rows.append(row)
        print(f"DF {kind:6s} {f / 1e9:4.1f} GHz  {row}", flush=True)
    return rows


# ---- 3. the analog beam with delay-line states ------------------------------------------------------
def analog_beam(f, theta, phi=0.0):
    """Realised gain loss (dB) against exact phases, and the peak sidelobe (dB), for a beam the
    2-bit delay lines steer toward theta at frequency f (identical cos^1.2 elements)."""
    k = 2 * np.pi * f / B.C0
    u0, v0 = B.direction(theta, phi)
    s = phase_states(u0, v0, f)
    w = np.exp(-1j * state_lag(s, f))
    th = np.linspace(-89.5, 89.5, 717)
    t, p = np.radians(th), np.radians(phi)
    u, v = np.sin(t) * np.cos(p), np.sin(t) * np.sin(p)
    af = np.abs(np.exp(1j * k * (u[:, None] * B.X + v[:, None] * B.Y)) @ w) ** 2 * np.cos(t) ** B.ELEM_Q
    peak_ideal = 16 ** 2 * np.cos(np.radians(theta)) ** B.ELEM_Q
    a0 = abs(np.sum(w * np.exp(1j * k * (B.X * u0 + B.Y * v0)))) ** 2 * np.cos(np.radians(theta)) ** B.ELEM_Q
    i = int(np.argmax(af))
    lo, hi = i, i
    while lo > 0 and af[lo - 1] < af[lo]:
        lo -= 1
    while hi < len(af) - 1 and af[hi + 1] < af[hi]:
        hi += 1
    rest = np.concatenate([af[:lo], af[hi + 1:]])
    return dict(loss_db=float(db(peak_ideal / a0)), sidelobe_db=float(db(rest.max() / af[i])) if len(rest) else -99.0,
                pointing_error_deg=float(th[i] - theta))


# ---- 4. accuracy limits: the fixed 4x4 against a resized circle --------------------------------------
def crlb_deg(positions, f, theta_deg, phi_deg, snr_db=10.0, samples=1000):
    """Cramer-Rao bound (one source, unknown amplitude, identical isotropic elements) on the angle
    theta in the plane phi, degrees: 1 / sqrt(2 K SNR sum (d psi_n / d theta - mean)^2)."""
    k = 2 * np.pi * f / B.C0
    x, y = positions
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    dpsi = k * (x * np.cos(p) + y * np.sin(p)) * np.cos(t)
    return float(np.degrees(1 / np.sqrt(2 * samples * 10 ** (snr_db / 10) * np.sum((dpsi - dpsi.mean()) ** 2))))


def uca(n, radius):
    b = 2 * np.pi * np.arange(n) / n
    return radius * np.cos(b), radius * np.sin(b)


def unambiguous_up_to(positions):
    """Highest frequency at which the closest pair of elements is still half a wavelength apart or less."""
    x, y = positions
    d = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    return float(B.C0 / (2 * d[d > 0].min()))


# ---- 5. other limits: grating lobes, the Chu bound, the Wilkinson combiners ---------------------------
def max_scan_without_grating_lobes(f, d=B.PITCH):
    s = B.C0 / (f * d) - 1
    return 90.0 if s >= 1 else float(np.degrees(np.arcsin(max(s, 0.0))))


def chu_bandwidth(f, radius):
    """Largest VSWR-2 fractional bandwidth of a lossless, single-tuned antenna inside a sphere of this
    radius (Chu / McLean Q = 1/(ka)^3 + 1/ka; BW = 1/(sqrt(2) Q))."""
    ka = 2 * np.pi * f / B.C0 * radius
    return float(1 / (np.sqrt(2) * (1 / ka ** 3 + 1 / ka)))


def wilkinson_tree(f, stages=2):
    """Lossless 4-way tree of quarter-wave (at F0) Wilkinson dividers: in-phase combining efficiency,
    input match and output isolation at f, by nodal analysis."""
    th = np.pi / 2 * f / F0
    z = np.sqrt(2) * Z0
    nodes = 1 + 2 + 4
    Y = np.zeros((nodes, nodes), complex)

    def line(i, j):
        yc, ys = -1j / (z * np.tan(th)), 1j / (z * np.sin(th))
        Y[i, i] += yc
        Y[j, j] += yc
        Y[i, j] += ys
        Y[j, i] += ys

    def res(i, j, r=2 * Z0):
        Y[i, i] += 1 / r
        Y[j, j] += 1 / r
        Y[i, j] -= 1 / r
        Y[j, i] -= 1 / r

    for i, j in ((0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)):
        line(i, j)
    for i, j in ((1, 2), (3, 4), (5, 6)):
        res(i, j)
    ports, inner = [0, 3, 4, 5, 6], [1, 2]
    Yp = Y[np.ix_(ports, ports)] - Y[np.ix_(ports, inner)] @ np.linalg.solve(Y[np.ix_(inner, inner)], Y[np.ix_(inner, ports)])
    S = np.linalg.solve(np.eye(5) + Z0 * Yp, np.eye(5) - Z0 * Yp)
    a = np.array([0, 0.5, 0.5, 0.5, 0.5])                   # equal in-phase waves into the four outputs
    return dict(combining_efficiency=float(abs(S[0] @ a) ** 2 / np.sum(abs(a) ** 2)),
                input_match_db=float(20 * np.log10(abs(S[0, 0]))), output_isolation_db=float(20 * np.log10(abs(S[1, 2]))))


# ---- 6. phase control for a 0.4-6 GHz board: which scheme serves DF and beams ----------------------
PHASE_OPTIONS = {   # available phase lags (radians) of one element's states at frequency f
    "delay lines, 90/180 deg at 5.8 GHz (today)": lambda f: state_lag(np.arange(4), f),
    "exact 0/180 switch (wideband balun)": lambda f: np.array([0.0, np.pi]),
    "exact 0/180 switch + 90-deg-at-5.8-GHz delay line": lambda f: np.array([0.0, np.pi, np.pi / 2 * f / F0, np.pi + np.pi / 2 * f / F0]),
    "3-bit phase shifter (45 deg steps at every f)": lambda f: np.arange(8) * np.pi / 4,
    "vector modulator (6 bits)": lambda f: np.arange(64) * np.pi / 32,
}


def quantised_beam(lags, positions, f, theta, phi=0.0):
    """Beam toward theta with each element on its nearest available phase (common phase searched):
    (gain loss against exact phases dB, highest sidelobe dB below the beam), cos^1.2 elements.
    Real (0/180-only) weights always give a mirror beam at -theta: a 0 dB 'sidelobe'."""
    k = 2 * np.pi * f / B.C0
    x, y = positions
    u0, v0 = B.direction(theta, phi)
    want = k * (x * u0 + y * v0)
    best, best_g = None, -1.0
    for c in np.linspace(0, 2 * np.pi, 96, endpoint=False):
        s = np.argmin(np.abs(np.angle(np.exp(1j * (want[:, None] + c - lags[None, :])))), axis=1)
        g = abs(np.sum(np.exp(1j * (want - lags[s]))))
        if g > best_g:
            best, best_g = lags[s], g
    w = np.exp(-1j * best)
    th = np.radians(np.linspace(-89.5, 89.5, 717))
    p = np.radians(phi)
    af = np.abs(np.exp(1j * k * ((np.sin(th) * np.cos(p))[:, None] * x + (np.sin(th) * np.sin(p))[:, None] * y)) @ w) ** 2
    af = af * np.cos(th) ** B.ELEM_Q
    i = int(np.argmax(af))
    lo, hi = i, i
    while lo > 0 and af[lo - 1] < af[lo]:
        lo -= 1
    while hi < len(af) - 1 and af[hi + 1] < af[hi]:
        hi += 1
    rest = np.concatenate([af[:lo], af[hi + 1:]])
    return float(db(len(x) ** 2 / best_g ** 2)), float(db(rest.max() / af[i])) if len(rest) else -99.0


def toggle_penalty_db(lags):
    """Decode noise penalty of the four-radio snapshot when the 'inversion' is the state closest to 180 degrees."""
    t = np.exp(-1j * lags[np.argmin(np.abs(np.angle(np.exp(1j * (lags - np.pi)))))])
    M = np.where(H4 < 0, t, 1)
    return float(db(np.trace(np.linalg.inv(M.conj().T @ M)).real))


def phase_control_options(freqs, scans=(0, 10, 20, 30, 40, 45)):
    """Each scheme on the built-in 4 x 4 panel (26 mm, so its beams lose meaning at low f) and on external
    antennas spaced half a wavelength at every band (4 x 4 grid)."""
    out = {}
    for name, lags_at in PHASE_OPTIONS.items():
        rows = []
        for f in freqs:
            lags = lags_at(f)
            ext = (B.X / B.PITCH * B.C0 / f / 2, B.Y / B.PITCH * B.C0 / f / 2)
            panel = [quantised_beam(lags, (B.X, B.Y), f, th) for th in scans]
            extn = [quantised_beam(lags, ext, f, th) for th in scans]
            rows.append(dict(f_GHz=round(f / 1e9, 3), decode_penalty_db=round(toggle_penalty_db(lags), 2),
                             panel_mean_loss_db=round(float(np.mean([r[0] for r in panel])), 2),
                             panel_worst_sidelobe_db=round(float(max(r[1] for r in panel)), 1),
                             external_mean_loss_db=round(float(np.mean([r[0] for r in extn])), 2),
                             external_worst_sidelobe_db=round(float(max(r[1] for r in extn)), 1)))
        out[name] = rows
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {"note": "simulator RF parts are frequency-flat: off 5.8 GHz this is the architecture with wideband parts"}

    f_fine = np.linspace(0.3e9, 11.3e9, 221)
    res["decode_noise_gain_db"] = [dict(f_GHz=round(f / 1e9, 2), toggle_deg=round(float(np.degrees(np.pi * f / F0)), 1),
                                        **{k: round(float(db(v)), 2) for k, v in decode_noise_gain(f).items()}) for f in f_fine]

    directions = [(float(t), float(p)) for t in (10, 25, 40) for p in range(0, 360, 30)]
    freqs = [0.9e9, 1.3e9, 1.8e9, 2.4e9, 3.0e9, 3.5e9, 4.5e9, 5.2e9, 5.8e9]
    res["df_vs_frequency"] = {kind: df_vs_frequency(freqs, kind, directions, compare_old=(1.8e9, 2.4e9, 3.5e9, 4.5e9))
                              for kind in ("quad", "hackrf")}

    beam_f = np.linspace(1.0e9, 6.5e9, 45)
    res["analog_beam"] = {str(th): [dict(f_GHz=round(f / 1e9, 3), **analog_beam(f, th)) for f in beam_f] for th in (0, 15, 30, 45)}

    grid = (B.X, B.Y)
    kraken_f = np.geomspace(0.1e9, 6e9, 60)
    res["accuracy_crlb_deg"] = [dict(
        f_GHz=round(float(f / 1e9), 3),
        array16_4x4=crlb_deg(grid, f, 20, 0),
        uca5_resized_spacing_0p5=crlb_deg(uca(5, 0.5 * B.C0 / f / (2 * np.sin(np.pi / 5))), f, 20, 0),
        uca5_fixed_for_1GHz=crlb_deg(uca(5, 0.5 * B.C0 / 1e9 / (2 * np.sin(np.pi / 5))), f, 20, 0),
        array16_beamwidth_deg=float(min(180.0, np.degrees(0.886 * B.C0 / f / (4 * B.PITCH)))) if B.C0 / f / (4 * B.PITCH) < 1.1 else 180.0,
    ) for f in kraken_f]
    res["unambiguous_up_to_GHz"] = dict(array16_4x4=round(unambiguous_up_to(grid) / 1e9, 2),
                                        uca5_fixed_for_1GHz=round(unambiguous_up_to(uca(5, 0.5 * B.C0 / 1e9 / (2 * np.sin(np.pi / 5)))) / 1e9, 2))

    lim_f = np.linspace(0.3e9, 8e9, 155)
    cell_radius = B.PITCH / np.sqrt(2)
    res["limits"] = [dict(f_GHz=round(f / 1e9, 3), max_scan_deg=max_scan_without_grating_lobes(f),
                          chu_bw_percent_26mm_cell=100 * chu_bandwidth(f, cell_radius),
                          chu_bw_percent_whole_array=100 * chu_bandwidth(f, 2 * B.PITCH * np.sqrt(2)),
                          **{"wilkinson_" + k: v for k, v in wilkinson_tree(f).items()}) for f in lim_f]
    res["phase_control_options"] = phase_control_options([0.4e9, 0.6e9, 0.9e9, 1.2e9, 1.5e9, 2.0e9, 2.4e9, 3.0e9, 3.5e9,
                                                          4.0e9, 4.5e9, 5.2e9, 5.8e9])
    patch = A.design(A.Substrate("FR-4, ground on L3", 4.4, 0.018, 1.275e-3, rough=1.8), probe_r=0.15e-3)
    iso = A.PatchArray(patch, [0.0], [0.0])
    res["patch_s11_db"] = [dict(f_GHz=round(f / 1e9, 3), s11_db=float(20 * np.log10(abs(iso.solve(f)["S"][0, 0])))) for f in lim_f]
    (Path("experiments") / "wideband_results.json").write_text(json.dumps(res, indent=1))

    # ---- figures ---------------------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 4))
    ng = res["decode_noise_gain_db"]
    fx = [r["f_GHz"] for r in ng]
    ax.plot(fx, [r["hackrf_16"] for r in ng], label="HackRF mode (16 patterns)")
    ax.plot(fx, [r["quad_4"] for r in ng], "--", label="four radios (4 patterns per quadrant)")
    ax.axvline(5.8, color="0.5", lw=0.8, ls=":")
    ax.set(xlabel="GHz", ylabel="noise after decoding, dB (0 = exact 180-degree flips)", ylim=(-1, 30),
           title="The 180-degree delay line away from 5.8 GHz: decode noise penalty")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "decode_noise.png", dpi=130)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for kind, style in (("quad", "o-"), ("hackrf", "s-")):
        rows = res["df_vs_frequency"][kind]
        ax.plot([r["f_GHz"] for r in rows], [r["frequency_aware"]["p95"] for r in rows], style, label=f"{kind}: frequency-aware decode, 95th pct")
        old = [r for r in rows if "decode_assuming_180" in r]
        ax.plot([r["f_GHz"] for r in old], [r["decode_assuming_180"]["p95"] for r in old], style.replace("-", ":"), alpha=0.6,
                label=f"{kind}: decode assuming exact 180, 95th pct")
    ax.set(xlabel="GHz", ylabel="direction error, degrees", yscale="log",
           title="Direction finding against frequency (wideband RF assumed, calibrated per frequency)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "df_vs_frequency.png", dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for th, rows in res["analog_beam"].items():
        fx = [r["f_GHz"] for r in rows]
        axs[0].plot(fx, [r["loss_db"] for r in rows], label=f"scan {th} deg")
        axs[1].plot(fx, [r["sidelobe_db"] for r in rows], label=f"scan {th} deg")
    axs[0].set(xlabel="GHz", ylabel="gain lost to the phase states, dB", title="Analog beam with the 2-bit delay lines")
    axs[1].set(xlabel="GHz", ylabel="highest sidelobe, dB below the beam", title="... and its sidelobes", ylim=(-20, 0))
    for ax in axs:
        ax.axvline(5.8, color="0.5", lw=0.8, ls=":")
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "analog_beam.png", dpi=130)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cr = res["accuracy_crlb_deg"]
    fx = [r["f_GHz"] for r in cr]
    ax.loglog(fx, [r["array16_4x4"] for r in cr], label="array16, fixed 4 x 4 at 26 mm")
    ax.loglog(fx, [r["uca5_resized_spacing_0p5"] for r in cr], label="5-element circle resized per band (KrakenSDR, spacing 0.5 wavelength)")
    ax.loglog(fx, [r["uca5_fixed_for_1GHz"] for r in cr], ":", label="5-element circle left at its 1 GHz size")
    ax.axvline(res["unambiguous_up_to_GHz"]["uca5_fixed_for_1GHz"], color="tab:green", lw=0.8, ls="--")
    ax.axvline(res["unambiguous_up_to_GHz"]["array16_4x4"], color="tab:blue", lw=0.8, ls="--")
    ax.set(xlabel="GHz", ylabel="best possible angle error, degrees (CRLB)",
           title="Noise-limited accuracy: per-element SNR 10 dB, 1000 samples\n(dashed: highest unambiguous frequency)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "accuracy_vs_aperture.png", dpi=130)

    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    lm = res["limits"]
    fx = [r["f_GHz"] for r in lm]
    axs[0, 0].plot(fx, [r["max_scan_deg"] for r in lm])
    axs[0, 0].set(title="26 mm pitch: widest scan without grating lobes", ylabel="degrees", xlabel="GHz")
    axs[0, 1].semilogy(fx, [r["chu_bw_percent_26mm_cell"] for r in lm], label="one element in a 26 mm cell")
    axs[0, 1].semilogy(fx, [r["chu_bw_percent_whole_array"] for r in lm], label="the whole 104 mm array")
    axs[0, 1].set(title="Chu limit: widest single-tuned bandwidth (VSWR 2)", ylabel="percent", xlabel="GHz", ylim=(0.1, 300))
    axs[0, 1].legend(fontsize=7)
    axs[1, 0].plot(fx, [10 * np.log10(r["wilkinson_combining_efficiency"]) for r in lm], label="in-phase combining loss")
    axs[1, 0].plot(fx, [r["wilkinson_output_isolation_db"] for r in lm], label="isolation between neighbouring inputs")
    axs[1, 0].set(title="4-way Wilkinson tree tuned to 5.8 GHz", ylabel="dB", xlabel="GHz", ylim=(-40, 1))
    axs[1, 0].legend(fontsize=7)
    axs[1, 1].plot(fx, [r["s11_db"] for r in res["patch_s11_db"]])
    axs[1, 1].axhline(-10, color="0.5", lw=0.8, ls=":")
    axs[1, 1].set(title="The FR-4 patch's match", ylabel="|S11|, dB", xlabel="GHz", ylim=(-30, 0))
    for ax in axs.flat:
        ax.axvline(5.8, color="0.5", lw=0.8, ls=":")
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "band_limits.png", dpi=130)

    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    for name, rows in res["phase_control_options"].items():
        fx = [r["f_GHz"] for r in rows]
        axs[0].plot(fx, [r["decode_penalty_db"] for r in rows], "o-", ms=3, label=name)
        axs[1].plot(fx, [r["external_mean_loss_db"] for r in rows], "o-", ms=3)
        axs[2].plot(fx, [r["external_worst_sidelobe_db"] for r in rows], "o-", ms=3)
    axs[0].set(title="Direction finding: decode noise penalty", ylabel="dB", xlabel="GHz")
    axs[1].set(title="Beam: gain lost to phase quantisation\n(antennas half a wavelength apart)", ylabel="dB", xlabel="GHz")
    axs[2].set(title="Beam: worst sidelobe, scans 0-45 deg", ylabel="dB below the beam", xlabel="GHz", ylim=(-14, 0.5))
    for ax in axs:
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(OUT / "phase_control.png", dpi=130)
    print("wrote experiments/wideband_results.json and", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))


if __name__ == "__main__":
    main()
