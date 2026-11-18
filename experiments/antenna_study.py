"""Antenna study: element options, coupling, scanning and direction finding with real patches.

    python experiments/antenna_study.py

Writes experiments/antenna_results.json and the figures in experiments/antenna/. The model and
its limits are described in array16/antenna.py; the write-up is experiments/antenna_study.md.
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, ".")
from array16 import Array, antenna as A, board as B, dsp  # noqa: E402
from array16.control import phase_states  # noqa: E402
from array16.sim import BoardModel, SimBackend, Source  # noqa: E402

OUT = Path("experiments/antenna")
F0, K0 = B.F0, 2 * np.pi * B.F0 / B.C0
ISM = (5.725e9, 5.875e9)
CENTRE, E_NB, H_NB, DIAG = 5, 6, 9, 10            # an inner element and its E-plane, H-plane and diagonal neighbours
SHOWN = {"centre (5)": 5, "edge, E-plane side (4)": 4, "edge, H-plane side (1)": 1, "corner (0)": 0}
PLANES = {"E": (0, 180), "H": (90, 270), "D": (45, 225)}   # both halves of each principal plane

OPTIONS = {   # key: (substrate, design keywords, eps_r tolerance for the "built board differs" case)
    "fr4_l2": (A.Substrate("FR-4 main board, ground on L2 (0.21 mm)", 4.4, 0.020, 0.21e-3, rough=1.8), dict(probe_r=0.15e-3), 0.2),
    "fr4_l3": (A.Substrate("FR-4 main board, ground on L3 (1.28 mm)", 4.4, 0.018, 1.275e-3, rough=1.8), dict(probe_r=0.15e-3), 0.2),
    "fr4_sep": (A.Substrate("FR-4 antenna board, 1.55 mm", 4.4, 0.018, 1.55e-3, rough=1.8), dict(probe_r=0.25e-3), 0.2),
    "ro4003_081": (A.Substrate("RO4003C antenna board, 0.813 mm", 3.55, 0.0027, 0.813e-3, rough=1.3), dict(probe_r=0.25e-3), 0.05),
    "ro4003_152": (A.Substrate("RO4003C antenna board, 1.524 mm", 3.55, 0.0027, 1.524e-3, rough=1.3), dict(probe_r=0.25e-3), 0.05),
    "ptfe_15": (A.Substrate("PTFE (F4B) antenna board, 1.5 mm", 2.65, 0.001, 1.5e-3, rough=1.2), dict(probe_r=0.25e-3), 0.05),
    "air_3": (A.Substrate("Sheet-metal patches 3 mm above the ground", 1.0, 0.0, 3e-3), dict(probe_r=0.5e-3, probe_comp=True), 0.0),
    "air_5": (A.Substrate("Sheet-metal patches 5 mm above the ground", 1.0, 0.0, 5e-3), dict(probe_r=0.5e-3, probe_comp=True), 0.0),
}
DF_OPTIONS = ["fr4_l3", "ro4003_152", "ptfe_15"]
# The slot model does not hold for the air options at 26 mm (facing edges 3.7 and 5.9 mm apart); these are
# the review's moment-method results (rooftop MPIE, exact image Green's function in air), 2026-11-15.
AIR_MOM = {
    "air_5": dict(e_plane_db=-9.9, e_plane_edge_db=-9.5, h_plane_db=-18.2, diagonal_db=-21.7, max_active_gamma_to_45deg=0.67,
                  df_boresight_cal_deg=dict(median=2.47, p95=5.71, worst=6.80), e_scan_loss_db_30_45_60=[0.06, 1.58, 6.56]),
    "air_3": dict(e_plane_db=-9.5, e_plane_edge_db=-8.7, max_active_gamma_to_45deg=0.72),
}
COLORS = dict(fr4_l2="#9467bd", fr4_l3="#d62728", fr4_sep="#ff7f0e", ro4003_081="#8c564b", ro4003_152="#1f77b4",
              ptfe_15="#17becf", air_3="#7f7f7f", air_5="#2ca02c", ideal="#000000")


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-30))


def cut(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.sin(t) * np.cos(p), np.sin(t) * np.sin(p)


def weights(arr, u0, v0, how):
    """Phase-only weights steering toward (u0, v0).
    manifold: conjugate of the true embedded patterns (the best a phase-only beam can do);
    boresight: ideal plane-wave phases after a boresight calibration, phases exact;
    2bit: the same rounded to the 2-bit states: HackRF mode, and the four-radio receive beam of
          Array.beam(word) without `toward`;
    2bit_quad: 2-bit inside each quadrant, each quadrant channel phased toward the beam: four-radio
          transmit (point(tx=True)), monopulse, and Array.beam(word, toward=...)."""
    e0 = arr.eep(np.array(0.0), np.array(0.0), F0)
    if how == "manifold":
        e = arr.eep(np.array(u0), np.array(v0), F0)
        return np.conj(e) / np.abs(e)
    if how == "boresight":
        return np.exp(-1j * K0 * (B.X * u0 + B.Y * v0)) * np.conj(e0) / np.abs(e0)
    w = np.exp(-0.5j * np.pi * phase_states(u0, v0, F0, np.angle(e0)))
    if how == "2bit_quad":
        ideal = np.exp(1j * K0 * (B.X * u0 + B.Y * v0)) * e0 / np.abs(e0)       # what the calibrated board expects
        for q in B.QUAD_ELEMS:
            w[q] *= np.exp(-1j * np.angle(np.sum(w[q] * ideal[q])))
    return w


def realised_gain(arr, w, u, v):
    return np.abs(arr.eep(np.asarray(u), np.asarray(v), F0) @ w) ** 2 / np.sum(np.abs(w) ** 2)


def scan_curves(arr, phi, thetas):
    out = {}
    for how in ("manifold", "boresight", "2bit", "2bit_quad"):
        g, gam = [], []
        for th in thetas:
            u0, v0 = B.direction(th, phi)
            w = weights(arr, u0, v0, how)
            g.append(float(db(realised_gain(arr, w, u0, v0))))
            gam.append(float(np.max(np.abs(arr.solve(F0)["S"] @ w / w))))
        out[how] = dict(gain_dbi=g, max_active_gamma=gam)
    return out


def peak_sidelobe(arr, w, phi, th0):
    """Highest lobe outside the main beam on the scan cut, dB below the beam."""
    th = np.linspace(-89.5, 89.5, 1433)
    p = realised_gain(arr, w, *cut(th, phi))
    i = int(np.argmin(np.abs(th - th0)))
    while 0 < i < len(p) - 1 and max(p[i - 1], p[i + 1]) > p[i]:       # climb to the beam's peak first
        i += 1 if p[i + 1] > p[i - 1] else -1
    lo, hi = i, i
    while lo > 0 and p[lo - 1] < p[lo]:
        lo -= 1
    while hi < len(p) - 1 and p[hi + 1] < p[hi]:
        hi += 1
    rest = np.concatenate([p[:lo], p[hi + 1:]])
    return float(db(rest.max() / p[i])) if len(rest) else float("-inf")


def df_errors(built, directions, design=None, seeds=(11, 12)):
    """Direction finding in four-radio mode on a board whose antenna is `built`, calibrated on a
    boresight source as calibrate() does. The same snapshots are also located against
    - the exact response of `built` (a bound: the model that made the snapshots), and
    - the simulated response of `design`, the nominal model, when the built board differs from it."""
    out = {"boresight_cal": [], "exact_manifold": [], "design_manifold": []}
    steer = {}
    if built is not None:
        e0 = built.eep(np.array(0.0), np.array(0.0), F0)
        steer["exact_manifold"] = lambda uu, vv: built.eep(uu, vv, F0) / e0
    if design is not None:
        d0 = design.eep(np.array(0.0), np.array(0.0), F0)
        steer["design_manifold"] = lambda uu, vv: design.eep(uu, vv, F0) / d0
    for seed in seeds:
        b = SimBackend("quad", seed=seed, board=BoardModel(seed + 1, antenna=built))
        a = Array(b)
        b.sources = [Source(0, 0, -75)]
        a.calibrate()
        for d in directions:
            b.sources = [Source(*d, power_dbm=-75)]
            z = a.snapshot()
            out["boresight_cal"].append(dsp.angle_between(B.angles(*dsp.locate(z, F0)), d))
            for name, s in steer.items():
                out[name].append(dsp.angle_between(B.angles(*dsp.locate(z, F0, s)), d))
    st = lambda e: dict(median=round(float(np.median(e)), 2), p95=round(float(np.percentile(e, 95)), 2), worst=round(float(np.max(e)), 2))
    return {k: st(v) for k, v in out.items() if v}, out["boresight_cal"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {"model": "array16/antenna.py: transmission-line patches coupled through space; see its docstring for accuracy"}
    arrays, curves, elements = {}, {}, {}
    thetas = list(range(0, 61, 5))
    upto45 = slice(0, thetas.index(45) + 1)

    # ---- every option: element, coupling, broadside array, scanning ------------------------------
    for key, (sub, kw, _) in OPTIONS.items():
        p = A.design(sub, **kw)
        c = A.characterise(p)
        arr = arrays[key] = A.PatchArray(p)
        S = arr.solve(F0)["S"]
        e0 = arr.eep(np.array(0.0), np.array(0.0), F0)
        c.update(substrate=sub.name, eps_r=sub.eps_r, tan_d=sub.tan_d, h_mm=sub.h * 1e3,
                 feed_from_edge_mm=(p.L / 2 - abs(p.feed)) * 1e3, coupling_reliable=arr.coupling_reliable)
        if arr.coupling_reliable:
            C, resid = A.fit_coupling(arr)
            Cn = C / np.mean(np.diag(C))
            c.update(
                coupling_db=dict(e_plane=float(20 * np.log10(abs(S[CENTRE, E_NB]))), h_plane=float(20 * np.log10(abs(S[CENTRE, H_NB]))),
                                 diagonal=float(20 * np.log10(abs(S[CENTRE, DIAG])))),
                receive_coupling_C=dict(e_plane_db=float(20 * np.log10(abs(Cn[CENTRE, E_NB]))), e_plane_deg=float(np.degrees(np.angle(Cn[CENTRE, E_NB]))),
                                        h_plane_db=float(20 * np.log10(abs(Cn[CENTRE, H_NB]))), h_plane_deg=float(np.degrees(np.angle(Cn[CENTRE, H_NB]))),
                                        fit_residual=resid),
                embedded_gain_dbi={n: float(db(abs(e0[i]) ** 2)) for n, i in SHOWN.items()},
                array_gain_broadside_dbi=float(db(abs(e0.sum()) ** 2 / 16)),
            )
            sc = {ph: scan_curves(arr, ph, thetas) for pair in PLANES.values() for ph in pair}
            curves[key] = sc
            worst = {pl: {how: np.minimum(sc[a][how]["gain_dbi"], sc[b][how]["gain_dbi"]) for how in sc[a]} for pl, (a, b) in PLANES.items()}
            g0 = worst["E"]["manifold"][0]
            c["scan_loss_db_worse_side"] = {f"{pl}_{th}": round(float(g0 - worst[pl]["manifold"][thetas.index(th)]), 2)
                                           for pl in PLANES for th in (30, 45, 60)}
            for how in ("2bit", "2bit_quad"):
                c[f"{how}_extra_loss_db_0_to_45deg"] = {pl: dict(mean=round(float(np.mean(worst[pl]["manifold"][upto45] - worst[pl][how][upto45])), 2),
                                                                 worst=round(float(np.max(worst[pl]["manifold"][upto45] - worst[pl][how][upto45])), 2))
                                                        for pl in PLANES}
            c["max_active_gamma_to_45deg"] = round(max(max(sc[ph]["manifold"]["max_active_gamma"][upto45]) for pair in PLANES.values() for ph in pair), 3)
            c["peak_sidelobe_db_worse_side"] = {f"{pl}_{th}_{how}": round(max(peak_sidelobe(arr, weights(arr, *B.direction(th, ph), how), ph, th)
                                                                              for ph in PLANES[pl]), 1)
                                                for pl in ("E", "H") for th in (0, 20, 40) for how in ("manifold", "2bit", "2bit_quad")}
        else:
            c["moment_method_reference"] = AIR_MOM[key]
        elements[key] = {k: (round(v, 3) if isinstance(v, float) else v) for k, v in c.items()}
        print(f"{key:11s} BW {c['bw_10dB_MHz']:4.0f} MHz  eff {c['efficiency']:.2f}  reliable {arr.coupling_reliable}", flush=True)
    res["options"] = elements

    # ---- direction finding -----------------------------------------------------------------------------
    directions = [(float(t), float(p)) for t in (5, 15, 25, 35, 45, 55) for p in range(0, 360, 30)]
    df, per_theta = {}, {}
    st, per_theta["ideal"] = df_errors(None, directions)
    df["ideal elements (simulator default)"] = st
    for key in DF_OPTIONS:
        sub, kw, tol = OPTIONS[key]
        nominal = arrays[key]
        built = A.PatchArray(replace(nominal.p, sub=replace(sub, eps_r=sub.eps_r + tol)))      # same layout, laminate off by tol
        st_nom, per_theta[key] = df_errors(nominal, directions)
        st_off, _ = df_errors(built, directions, design=nominal)
        df[key] = dict(as_designed=st_nom, eps_r_off=dict(delta_eps_r=tol, **st_off))
        print(f"DF {key:11s} as designed {st_nom}\n   eps_r +{tol}: {st_off}", flush=True)
    res["direction_finding_deg"] = df

    # ---- how much coupling a boresight calibration tolerates: the FR-4 coupling matrix, scaled and turned ---
    C, _ = A.fit_coupling(arrays["fr4_l3"])
    Cn = C / np.mean(np.diag(C))
    Dg, O = np.diag(np.diag(Cn)), Cn - np.diag(np.diag(Cn))
    ref = abs(O[CENTRE, E_NB])
    sens = []
    for lvl in (-30, -25, -20, -17, -15, -12):
        runs = []
        for psi in (0, 90, 180, 270):
            model = A.CoupledIdeal(C=Dg + 10 ** (lvl / 20) / ref * np.exp(1j * np.radians(psi)) * O)
            runs.append(df_errors(model, directions, seeds=(11,))[0])
        row = dict(e_plane_coupling_db=lvl,
                   boresight_cal_median=[min(r["boresight_cal"]["median"] for r in runs), max(r["boresight_cal"]["median"] for r in runs)],
                   boresight_cal_p95=[min(r["boresight_cal"]["p95"] for r in runs), max(r["boresight_cal"]["p95"] for r in runs)],
                   exact_manifold_p95=[min(r["exact_manifold"]["p95"] for r in runs), max(r["exact_manifold"]["p95"] for r in runs)])
        sens.append(row)
        print(f"coupling {lvl:4d} dB: {row}", flush=True)
    res["df_vs_coupling_level"] = dict(note="FR-4 receive-coupling matrix (all neighbours) scaled to the level and its phase turned 0/90/180/270 deg; [min, max] over phase",
                                       rows=sens)
    (Path("experiments") / "antenna_results.json").write_text(json.dumps(res, indent=1))

    # ---- figures ---------------------------------------------------------------------------------------
    f = np.linspace(5.3e9, 6.3e9, 401)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for key, arr in arrays.items():
        iso = A.PatchArray(arr.p, [0.0], [0.0])
        ax.plot(f / 1e9, [20 * np.log10(abs(iso.solve(x)["S"][0, 0])) for x in f], color=COLORS[key], label=OPTIONS[key][0].name)
    ax.axvspan(ISM[0] / 1e9, ISM[1] / 1e9, color="0.85", zorder=0)
    ax.axhline(-10, color="0.4", lw=0.8, ls=":")
    ax.set(xlabel="GHz", ylabel="|S11| of one patch, dB", ylim=(-35, 0), title="Element match (ISM band shaded)")
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "element_match.png", dpi=130)

    th = np.linspace(-90, 90, 361)
    fig, axs = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for row, key in enumerate(("fr4_l3", "ptfe_15")):
        arr = arrays[key]
        iso = A.PatchArray(arr.p, [0.0], [0.0])
        for col, (plane, phi) in enumerate((("E-plane", 0), ("H-plane", 90))):
            ax = axs[row, col]
            u, v = cut(th, phi)
            e = arr.eep(u, v, F0)
            for name, i in SHOWN.items():
                ax.plot(th, db(np.abs(e[:, i]) ** 2), label=name)
            ax.plot(th, db(np.abs(iso.eep(u, v, F0)[:, 0]) ** 2), "k--", lw=1, label="same patch alone")
            ax.set(title=f"{OPTIONS[key][0].name}: {plane}", ylim=(-25, 10), ylabel="realised gain, dBi")
            ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=7)
    for ax in axs[1]:
        ax.set_xlabel("theta, degrees")
    fig.tight_layout()
    fig.savefig(OUT / "embedded_patterns.png", dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for col, key in enumerate(("fr4_l3", "ptfe_15")):
        arr, ax = arrays[key], axs[col]
        u, v = cut(th, 0)
        e, e0 = arr.eep(u, v, F0), arr.eep(np.array(0.0), np.array(0.0), F0)
        ideal = np.exp(1j * K0 * (u[:, None] * B.X + v[:, None] * B.Y))
        r = e / e0 / ideal
        err = np.degrees(np.angle(r * np.conj(np.exp(1j * np.angle(r.sum(axis=1))))[:, None]))     # less the phase common to all 16
        for name, i in SHOWN.items():
            ax.plot(th, err[:, i], label=name)
        ax.set(title=f"{OPTIONS[key][0].name}\nphase left after a boresight calibration (E-plane, common phase removed)",
               xlabel="theta, degrees", ylabel="degrees", ylim=(-40, 40))
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "calibration_residue.png", dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for col, plane in enumerate(("E", "H")):
        ax = axs[col]
        a, b = PLANES[plane]
        for key in ("fr4_l2", "fr4_l3", "ro4003_152", "ptfe_15"):
            for how, ls in (("manifold", "-"), ("2bit", "--")):
                ax.plot(thetas, np.minimum(curves[key][a][how]["gain_dbi"], curves[key][b][how]["gain_dbi"]), color=COLORS[key], ls=ls,
                        lw=1.6 if ls == "-" else 1, label=OPTIONS[key][0].name if how == "manifold" else None)
        ax.set(title=f"{plane}-plane scan, worse side", xlabel="scan angle, degrees")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("realised gain, dBi")
    axs[0].legend(fontsize=7)
    fig.suptitle("Array gain toward the beam (solid: exact phases; dashed: the hardware's 2-bit states)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "scan_gain.png", dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    ths = np.array([d[0] for d in directions] * 2)
    for key, errs in per_theta.items():
        label = "ideal elements" if key == "ideal" else OPTIONS[key][0].name
        axs[0].plot(sorted(errs), np.linspace(0, 1, len(errs)), color=COLORS[key], label=label)
        means = [np.mean(np.array(errs)[ths == t]) for t in (5, 15, 25, 35, 45, 55)]
        axs[1].plot((5, 15, 25, 35, 45, 55), means, "o-", color=COLORS[key])
    axs[0].set(xscale="log", xlabel="direction error, degrees", ylabel="fraction of directions",
               title="Direction finding, boresight calibration only")
    axs[1].set(xlabel="source theta, degrees", ylabel="mean error, degrees", title="Error grows off boresight")
    axs[0].legend(fontsize=7)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "df_error.png", dpi=130)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    lv = [s["e_plane_coupling_db"] for s in sens]
    ax.fill_between(lv, [s["boresight_cal_p95"][0] for s in sens], [s["boresight_cal_p95"][1] for s in sens], alpha=0.3,
                    label="boresight calibration, 95th pct (band: coupling phase)")
    ax.plot(lv, [np.mean(s["boresight_cal_median"]) for s in sens], "o-", label="boresight calibration, median")
    ax.fill_between(lv, [s["exact_manifold_p95"][0] for s in sens], [s["exact_manifold_p95"][1] for s in sens], alpha=0.3, color="tab:green",
                    label="exact array response, 95th pct (bound)")
    ax.set(xlabel="nearest-neighbour (E-plane) coupling, dB", ylabel="direction error, degrees", ylim=(0, None),
           title="How much coupling a boresight calibration tolerates")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "df_vs_coupling.png", dpi=130)
    print("wrote experiments/antenna_results.json and", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))


if __name__ == "__main__":
    main()
