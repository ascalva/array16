"""Virtual bring-up: run the bench procedure on many simulated boards and print the statistics.

    python experiments/virtual_bringup.py [trials]
"""
import json
import sys

import numpy as np

sys.path.insert(0, ".")
from array16 import Array, board as B, dsp                      # noqa: E402
from array16.board import ControlWord                           # noqa: E402
from array16.sim import BoardModel, SimBackend, Source          # noqa: E402

IDEAL = dict(amp_sigma_db=0, phase_sigma_deg=0, stage_sigma_deg=0, bypass_sigma_deg=0, port_sigma_deg=0)


def bench(kind, seed, ideal=False):
    b = SimBackend(kind, seed=seed, board=BoardModel(seed + 1, **IDEAL) if ideal else None)
    a = Array(b)
    if not ideal:
        b.sources = [Source(0, 0, -75)]
        a.calibrate()
        b.sources = []
    return b, a


def stats(err):
    return dict(median=round(float(np.median(err)), 2), p95=round(float(np.percentile(err, 95)), 2), worst=round(float(np.max(err)), 2))


def main(trials=40):
    rng = np.random.default_rng(2026)
    dirs = [(float(rng.uniform(5, 45)), float(rng.uniform(-180, 180))) for _ in range(trials)]
    out = {}

    # direction finding from an element snapshot
    cases = {"hackrf_tone": ("hackrf", "tone", False), "four_radio_tone": ("quad", "tone", False),
             "four_radio_modulated_blind": ("quad", "qpsk", True), "four_radio_modulated_naive": ("quad", "qpsk", False)}
    for name, (kind, sig, blind) in cases.items():
        err = []
        for i, d in enumerate(dirs):
            b, a = bench(kind, 100 + i)
            if kind == "quad":
                a.tune(B.F0)                                    # a retune before every measurement
            b.sources = [Source(*d, power_dbm=-75, kind=sig)]
            err.append(dsp.angle_between(a.locate(blind=blind), d))
        out["direction_error_deg_" + name] = stats(err)

    # single-burst direction from the four quadrants, source within 6 degrees of the steered beam
    err = []
    for i, d in enumerate(dirs):
        b, a = bench("quad", 300 + i)
        b.sources = [Source(d[0] + rng.uniform(-4, 4), d[1], -75, kind="qpsk")]
        err.append(dsp.angle_between(a.monopulse(*d), (b.sources[0].theta_deg, d[1])))
    out["direction_error_deg_single_burst"] = stats(err)

    # channel alignment after a retune
    res = []
    for i in range(trials):
        b, a = bench("quad", 500 + i)
        b.sources = [Source(0, 0, -70)]
        a.tune(B.F0)
        z = a.snapshot()
        res.append(np.degrees(np.std(np.angle(z * np.conj(z.mean())))))
    out["phase_error_after_retune_deg_rms"] = stats(res)

    # steered transmit beams on real boards: where the peak lands and how strong it is
    for name, kind, align in (("hackrf", "hackrf", False), ("four_radio_unaligned", "quad", False), ("four_radio_aligned", "quad", True)):
        point, eirp = [], []
        for i, d in enumerate(dirs):
            b, a = bench(kind, 700 + i)
            if align:
                a.calibrate_tx(lambda: 10 ** (b.eirp_dbm(0, 0) / 10))
            a.point(*d, tx=True)
            cut = np.array([b.eirp_dbm(th, d[1]) for th in np.arange(0, 70.5, 0.5)])
            point.append(abs(0.5 * int(np.argmax(cut)) - d[0]))
            eirp.append(b.eirp_dbm(*d))
        out["tx_pointing_error_deg_" + name] = stats(point)
        out["tx_eirp_dbm_steered_" + name] = dict(median=round(float(np.median(eirp)), 1), lowest=round(float(np.min(eirp)), 1))

    # link budget on an ideal board
    b, a = bench("quad", 1, ideal=True)
    a.calibrate_tx(lambda: 10 ** (b.eirp_dbm(0, 0) / 10))
    a.point(0, 0, tx=True)
    out["eirp_dbm_broadside_four_radio"] = round(b.eirp_dbm(0, 0), 2)
    b, a = bench("hackrf", 1, ideal=True)
    a.set_word(ControlWord.transmit(0, hackrf_mode=True))
    out["eirp_dbm_broadside_hackrf"] = round(b.eirp_dbm(0, 0), 2)
    for kind in ("quad", "hackrf"):
        b, a = bench(kind, 1, ideal=True)
        p_dbm = -105.0
        b.sources = [Source(0, 0, p_dbm, f_offset=20e3)]
        seg, _ = a._slots([ControlWord.receive(0, hackrf_mode=kind == "hackrf")] * 8)
        tone = np.exp(-2j * np.pi * 20e3 / b.fs * np.arange(seg.shape[-1]))
        amp = np.array([[np.mean(s[c] * tone) for c in range(seg.shape[2])] for s in seg[0]])
        sig = np.mean(np.abs(amp) ** 2, axis=0)
        noise = np.mean(np.abs(seg[0].transpose(1, 0, 2).reshape(seg.shape[2], -1)) ** 2, axis=1) - sig
        snr = seg.shape[2] * np.mean(sig / noise)
        kTB = 10 * np.log10(B.K_BOLTZ * B.T0 * b.fs) + 30
        out["noise_figure_db_" + kind] = round(float(p_dbm + B.ARRAY_GAIN_DBI - kTB - 10 * np.log10(snr)), 2)

    slot = 50e-6
    out["snapshot_ms"] = dict(hackrf=round(17 * slot * 1e3, 2), four_radio_steady=round(5 * slot * 1e3, 2), four_radio_any=round(6 * slot * 1e3, 2))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
