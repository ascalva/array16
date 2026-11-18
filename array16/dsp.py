"""Signal processing on captures: alignment, snapshot decoding, direction finding."""
from __future__ import annotations

import numpy as np

from . import board as B
from .control import BLIND_LINKS, H4, H16, toggle_factor

GUARD = 10            # samples dropped at each end of a slot (switching transients, marker jitter)


def cal_slot(y: np.ndarray, k0: int, slot: int, replica: np.ndarray, max_skew: int = 64):
    """Use the calibration slot that starts at sample k0 of board A.

    The probe sends a known pseudo-noise code. Each channel is correlated with
    the replica; the code phase found on board A (channels 0-1) against board B
    (channels 2-3) gives the whole-sample offset between the two streams, and
    the correlation peaks give the phase of every channel, with other signals
    in the slot suppressed by the correlation.
    Returns (offset of board B in samples, channel phasors relative to channel 0).
    """
    n, back = len(replica), 20                                  # board B's window sits 20 samples early
    lo, hi = k0 + GUARD, k0 + slot - GUARD
    rep = np.conj(np.fft.fft(replica))
    corr = lambda x: np.fft.ifft(np.fft.fft(x, n) * rep)        # corr[m] = sum x[k] conj(replica[k - m])
    c = [corr(y[0, lo:hi]), corr(y[1, lo:hi]), corr(y[2, lo - back:hi - back]), corr(y[3, lo - back:hi - back])]
    ma = int(np.argmax(np.abs(c[0]) ** 2 + np.abs(c[1]) ** 2))
    mb = int(np.argmax(np.abs(c[2]) ** 2 + np.abs(c[3]) ** 2))
    d = (ma - mb + back) % n
    d = d if d < max_skew else 0
    chan = np.array([c[0][ma], c[1][ma], c[2][mb], c[3][mb]])
    return d, chan / chan[0]


def shift_boards(y: np.ndarray, d: int) -> np.ndarray:
    """Apply the stream offset found by cal_slot so both boards share one time axis."""
    m = y.shape[1] - d
    return np.concatenate([y[:2, d:d + m], y[2:, :m]])


def find_marker(y: np.ndarray, slot: int, period: int) -> int:
    """Sample index where the blank marker slot starts (first occurrence)."""
    p = np.sum(np.abs(y) ** 2, axis=0)
    T = slot * period
    reps = len(p) // T
    folded = p[:reps * T].reshape(reps, T).mean(0)
    win = slot - 8
    c = np.concatenate([[0.0], np.cumsum(np.concatenate([folded, folded]))])
    means = (c[win:win + T] - c[:T]) / win
    return (int(np.argmin(means)) - 4) % T


def cut_slots(y: np.ndarray, marker: int, slot: int, period: int, reps: int = 2):
    """Slice `reps` repetitions into slots.

    Returns (seg, idx): seg[r, s] is the (channels, slot - 2*GUARD) block of data
    slot s in repetition r, and idx[r, s] the sample index of its first sample.
    """
    idx = np.array([[marker + (r * period + 1 + s) * slot + GUARD for s in range(period - 1)] for r in range(reps)])
    L = slot - 2 * GUARD
    if idx.max() + L > y.shape[1]:
        raise ValueError("capture too short for %d repetitions" % reps)
    seg = np.array([[y[:, i:i + L] for i in row] for row in idx])
    return seg, idx


def estimate_tone(seg: np.ndarray, idx: np.ndarray, fs: float) -> float:
    """Frequency of a steady tone in the slots, Hz from the tuned frequency.

    Three steps: a zero-padded spectrum, a half-slot phase slope, then the phase
    advance of each slot between two repetitions.
    """
    L = seg.shape[-1]
    spec = np.abs(np.fft.fft(seg, 4096, axis=-1)) ** 2
    f = np.fft.fftfreq(4096, 1 / fs)[np.argmax(spec.sum((0, 1, 2)))]
    k = np.arange(L)
    d = seg * np.exp(-2j * np.pi * f * k / fs)
    lag = L // 2
    f += np.angle(np.sum(d[..., lag:] * np.conj(d[..., :-lag]))) * fs / (2 * np.pi * lag)
    if seg.shape[0] > 1:
        m = slot_means(seg, idx, fs, f)
        f += np.angle(np.sum(m[1] * np.conj(m[0]))) * fs / (2 * np.pi * (idx[1, 0] - idx[0, 0]))
    return float(f)


def slot_means(seg: np.ndarray, idx: np.ndarray, fs: float, f: float) -> np.ndarray:
    """Mean of each slot after removing a tone at f. Phase is referred to sample 0 of the capture."""
    k = idx[..., None, None] + np.arange(seg.shape[-1])
    return np.mean(seg * np.exp(-2j * np.pi * f * k / fs), axis=-1)


def _patterns(H: np.ndarray, toggle: np.ndarray) -> np.ndarray:
    """Pattern matrix with each inverted entry replaced by that element's actual toggle factor."""
    return np.where(H < 0, toggle[None, :], 1)


def decode_coherent(g: np.ndarray, hackrf: bool, f: float = B.F0, base_state=0) -> np.ndarray:
    """Slot means g[slot, channel] -> one complex value per element (steady signals).

    At F0 an inverted element is exactly -1 and the Hadamard matrix inverts itself. At other
    frequencies the 180-degree line gives pi * f / F0 instead (control.toggle_factor) and the
    patterns are solved for; noise grows as that step moves away from 180 degrees.
    """
    if f == B.F0:
        if hackrf:
            return H16.T @ g[:, 0] / 16
        z = np.zeros(B.N_ELEM, complex)
        for q in range(4):
            z[B.QUAD_ELEMS[q]] = H4.T @ g[:, q] / 4
        return z
    t = toggle_factor(base_state, f)
    if hackrf:
        return np.linalg.solve(_patterns(H16, t), g[:, 0])
    z = np.zeros(B.N_ELEM, complex)
    for q in range(4):
        e = B.QUAD_ELEMS[q]
        z[e] = np.linalg.solve(_patterns(H4, t[e]), g[:, q])
    return z


def decode_blind(seg: np.ndarray, schedule, f: float = B.F0, base_state=0) -> np.ndarray:
    """Four-channel slots -> one complex value per element, for any signal.

    Each slot gives the four quadrant outputs up to an unknown common factor
    (the signal itself). A quadrant that holds pattern 0 across two slots ties
    their factors together; see control.BLIND_SCHEDULE. Off F0 the patterns are
    solved with the actual toggle factor, as in decode_coherent().
    """
    b = []
    for s in range(seg.shape[1]):
        yy = np.concatenate([seg[r, s] for r in range(seg.shape[0])], axis=-1)
        w, v = np.linalg.eigh(yy @ yy.conj().T / yy.shape[1])
        b.append(v[:, -1] * np.sqrt(max(w[-1] - w[:-1].mean(), 1e-30)))
    for s, q in enumerate(BLIND_LINKS):
        b[s + 1] = b[s + 1] * (b[s][q] / b[s + 1][q])
    z = np.zeros(B.N_ELEM, complex)
    t = None if f == B.F0 else toggle_factor(base_state, f)
    for q in range(4):
        g = np.array([np.mean([b[s][q] for s in range(len(schedule)) if schedule[s][q] == p]) for p in range(4)])
        e = B.QUAD_ELEMS[q]
        z[e] = H4.T @ g / 4 if t is None else np.linalg.solve(_patterns(H4, t[e]), g)
    return z


def locate(a: np.ndarray, f: float = B.F0, steer=None) -> tuple[float, float]:
    """Direction (u, v) of the strongest plane wave in an element snapshot.

    steer(u, v) -> (..., 16): the array's response toward (u, v) as the snapshot sees it, for
    example measured or simulated embedded element patterns divided by their boresight value
    (what calibrate() divides by). None: ideal plane-wave phases.
    """
    k = 2 * np.pi * f / B.C0

    def best(us, vs):
        uu, vv = np.meshgrid(us, vs)
        if steer is None:
            p = np.abs(np.exp(-1j * k * (uu[..., None] * B.X + vv[..., None] * B.Y)) @ a) ** 2
        else:
            s = steer(uu, vv)
            p = np.abs(np.conj(s) @ a) ** 2 / np.maximum(np.sum(np.abs(s) ** 2, axis=-1), 1e-30)
        p[uu ** 2 + vv ** 2 > 1] = 0
        i = np.unravel_index(np.argmax(p), p.shape)
        return float(uu[i]), float(vv[i])

    u, v = best(np.arange(-0.98, 0.99, 0.02), np.arange(-0.98, 0.99, 0.02))
    return best(u + np.arange(-0.03, 0.0301, 0.002), v + np.arange(-0.03, 0.0301, 0.002))


def angle_between(d1, d2) -> float:
    """Angle in degrees between two directions given as (theta_deg, phi_deg)."""
    def vec(d):
        u, v = B.direction(*d)
        return np.array([u, v, np.sqrt(max(0.0, 1 - u * u - v * v))])
    return float(np.degrees(np.arccos(np.clip(np.dot(vec(d1), vec(d2)), -1, 1))))
