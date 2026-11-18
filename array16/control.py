"""Beam and pattern generation: which control words to load for each job."""
from __future__ import annotations

import numpy as np

from .board import F0, C0, N_ELEM, QUAD_ELEMS, X, Y, ControlWord

H4 = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])
H16 = np.kron(H4, H4)

# schedule[slot][quadrant] = Hadamard pattern shown by that quadrant in that slot
COHERENT_SCHEDULE = [[k] * 4 for k in range(4)]
# Blind schedule: every quadrant shows all four patterns, and in each pair of
# neighbouring slots one quadrant holds pattern 0. That held pattern ties the
# unknown signal level of one slot to the next, so the signal may be modulated.
BLIND_SCHEDULE = [[0, 1, 1, 1], [0, 0, 2, 2], [1, 0, 0, 3], [2, 2, 0, 0], [3, 3, 3, 0]]
BLIND_LINKS = [0, 1, 2, 3]            # quadrant that holds pattern 0 across slots s and s + 1


def state_lag(state, f: float = F0) -> np.ndarray:
    """Phase lag of the 2-bit phase states at frequency f, radians.

    The phase shifter switches delay lines (true time delays): 90 and 180 degrees long at F0,
    so a state's lag grows in proportion to frequency and is a quarter turn per step only at F0.
    """
    return np.asarray(state) * (np.pi / 2) * (f / F0)


def phase_states(u: float, v: float, f: float = F0, offsets: np.ndarray | None = None) -> np.ndarray:
    """2-bit phase state (0..3) per element that points the beam at (u, v).

    offsets: measured phase of each element path at state 0, radians (from calibration).
    A common phase is free, so several values are tried and the best-combining one kept.
    Away from F0 the four states are not quarter turns apart (state_lag), and each element
    takes whichever state comes closest.
    """
    want = 2 * np.pi * f / C0 * (X * u + Y * v) + (0 if offsets is None else offsets)
    best, best_gain = None, -1.0
    if f == F0:
        for c in np.linspace(0, np.pi / 2, 32, endpoint=False):
            s = np.round((want + c) / (np.pi / 2)).astype(int)
            gain = abs(np.sum(np.exp(1j * (want - s * np.pi / 2))))
            if gain > best_gain:
                best, best_gain = s % 4, gain
        return best
    lags = state_lag(np.arange(4), f)
    for c in np.linspace(0, 2 * np.pi, 128, endpoint=False):
        s = np.argmin(np.abs(np.angle(np.exp(1j * (want[:, None] + c - lags[None, :])))), axis=1)
        gain = abs(np.sum(np.exp(1j * (want - lags[s]))))
        if gain > best_gain:
            best, best_gain = s, gain
    return best


def toggle_factor(base_state=0, f: float = F0) -> np.ndarray:
    """What switching each element's 180-degree line multiplies its signal by, at frequency f.

    Exactly -1 at F0; elsewhere a phase step of pi * f / F0 (a lag if the line is switched in).
    """
    s = np.broadcast_to(np.asarray(base_state, dtype=int), (N_ELEM,)) % 4
    return np.exp(-1j * (state_lag((s + 2) % 4, f) - state_lag(s, f)))


def snapshot_words(base_state=0, hackrf_mode=False, blind=False) -> tuple[list[ControlWord], list[list[int]] | None]:
    """Control words for an element snapshot.

    HackRF mode: 16 words, element n inverted in slot k when H16[k, n] < 0 (inverted means its
    180-degree line toggled: exactly 180 degrees at F0, see toggle_factor()).
    Four-radio mode: each quadrant steps through H4 on its own four elements,
    following COHERENT_SCHEDULE (steady signals) or BLIND_SCHEDULE (any signal).
    Returns (words, schedule); schedule is None in HackRF mode.
    """
    base = np.broadcast_to(np.asarray(base_state, dtype=int), (N_ELEM,))
    if hackrf_mode:
        return [ControlWord.receive((base + np.where(H16[k] < 0, 2, 0)) % 4, hackrf_mode=True) for k in range(16)], None
    schedule = BLIND_SCHEDULE if blind else COHERENT_SCHEDULE
    words = []
    for slot in schedule:
        s = base.copy()
        for q, pat in enumerate(slot):
            s[QUAD_ELEMS[q]] += np.where(H4[pat] < 0, 2, 0)
        words.append(ControlWord.receive(s % 4))
    return words, schedule


def scan_grid(max_theta_deg: float = 50.0, n: int = 8) -> list[tuple[float, float]]:
    """n x n grid of (u, v) directions inside max_theta, for scan mode."""
    lim = np.sin(np.radians(max_theta_deg))
    g = np.linspace(-lim, lim, n)
    return [(float(a), float(b)) for b in g for a in g if a * a + b * b <= lim * lim + 1e-9]
