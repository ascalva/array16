"""Board description: geometry, RF budget and the 56-bit control word.

The firmware model, the simulator, the DSP code and the schematic generator
all import their numbers from this file, so the hardware is described once.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

C0 = 299_792_458.0
F0 = 5.8e9                 # design centre frequency, Hz
PITCH = 26e-3              # element pitch, m
N_ELEM = 16
K_BOLTZ = 1.380649e-23
T0 = 290.0

_idx = np.arange(N_ELEM)
ROW, COL = _idx // 4, _idx % 4       # element e = 4*row + col; row 0 = top, col 0 = left, seen from the front
X = (COL - 1.5) * PITCH              # m, +x to the right
Y = (1.5 - ROW) * PITCH              # m, +y up
QUAD = (ROW // 2) * 2 + COL // 2     # 0 top-left, 1 top-right, 2 bottom-left, 3 bottom-right
QUAD_ELEMS = [np.flatnonzero(QUAD == q) for q in range(4)]

# ---- RF budget, dB unless noted (datasheets and the design doc) -------------
IL_SW = 0.46                         # BGS12WN6 insertion loss at 5.8 GHz
ISO_SW = 26.0                        # BGS12WN6 isolation
LNA_GAIN, LNA_NF = 13.0, 1.2         # SKY65404-31
LNA_OFF_GAIN = -25.0                 # SKY65404-31 disabled
LINE_DB_PER_CM = 0.32                # 50-ohm microstrip on the 0.21 mm prepreg
FEED_LOSS = 0.3                      # feed via and pad ahead of the first switch
WILK_EXCESS = 0.35                   # per Wilkinson stage, above the ideal 3.01 dB
LEN_PS_BASE_CM = 0.4                 # line common to every phase state
LEN_90_CM, LEN_180_CM = 0.724, 1.448 # extra line in the long path of each stage
LEN_QUAD_TREE_CM = 3.5               # line inside one quadrant's 4-way combiner
LEN_DIAG_CM = 3.7                    # quadrant port to the shared combiner
LEN_HACKRF_CM = 6.0                  # shared combiner to the HackRF socket
ARRAY_GAIN_DBI = 14.6                # realised broadside gain, all 16 elements
ELEM_Q = 1.2                         # embedded element power pattern ~ cos(theta)**q
RADIO_NF = {"quad": 6.0, "hackrf": 12.0}   # assumed, see the design doc

# ---- 56-bit control word = 7 x 74HC595 --------------------------------------
# Bit i of the word is output Q(i % 8) of shift register i // 8.
BIT_RX, BIT_P180, BIT_P90 = 0, 16, 32      # 16 bits each, indexed by element
BIT_MODE = 48                              # 1 = quadrants feed the shared combiner (HackRF mode)
BIT_QRX = 49                               # 4 bits: quadrant port receives (0 = transmit socket)
BIT_HRX = 53                               # HackRF port receives
BIT_CAL = 54                               # calibration tone routed to the centre probe
BIT_LNA2 = 55                              # second-stage LNAs enabled
WORD_BITS, WORD_BYTES = 56, 7


def _bits(n: int) -> np.ndarray:
    return np.zeros(n, dtype=bool)


@dataclass
class ControlWord:
    """One setting of every switch on the board. All zeros is the safe power-up state."""

    rx: np.ndarray = field(default_factory=lambda: _bits(N_ELEM))     # element receives through its LNA
    p180: np.ndarray = field(default_factory=lambda: _bits(N_ELEM))   # long path of the 180-degree stage
    p90: np.ndarray = field(default_factory=lambda: _bits(N_ELEM))    # long path of the 90-degree stage
    hackrf_mode: bool = False
    qrx: np.ndarray = field(default_factory=lambda: _bits(4))
    hrx: bool = False
    cal: bool = False
    lna2: bool = False

    @property
    def state(self) -> np.ndarray:
        """Phase lag per element in quarter turns (0..3)."""
        return 2 * self.p180.astype(int) + self.p90.astype(int)

    @classmethod
    def receive(cls, state=0, hackrf_mode=False, cal=False) -> "ControlWord":
        s = np.broadcast_to(np.asarray(state, dtype=int), (N_ELEM,)) % 4
        return cls(rx=np.ones(N_ELEM, bool), p180=s >= 2, p90=(s % 2) == 1, hackrf_mode=hackrf_mode,
                   qrx=np.full(4, not hackrf_mode), hrx=hackrf_mode, cal=cal, lna2=True)

    @classmethod
    def transmit(cls, state=0, hackrf_mode=False) -> "ControlWord":
        s = np.broadcast_to(np.asarray(state, dtype=int), (N_ELEM,)) % 4
        return cls(p180=s >= 2, p90=(s % 2) == 1, hackrf_mode=hackrf_mode)

    def to_int(self) -> int:
        v = 0
        for base, arr in ((BIT_RX, self.rx), (BIT_P180, self.p180), (BIT_P90, self.p90), (BIT_QRX, self.qrx)):
            for i, b in enumerate(arr):
                v |= int(bool(b)) << (base + i)
        for bit, b in ((BIT_MODE, self.hackrf_mode), (BIT_HRX, self.hrx), (BIT_CAL, self.cal), (BIT_LNA2, self.lna2)):
            v |= int(bool(b)) << bit
        return v

    @classmethod
    def from_int(cls, v: int) -> "ControlWord":
        get = lambda base, n: np.array([(v >> (base + i)) & 1 for i in range(n)], dtype=bool)
        one = lambda bit: bool((v >> bit) & 1)
        return cls(rx=get(BIT_RX, N_ELEM), p180=get(BIT_P180, N_ELEM), p90=get(BIT_P90, N_ELEM),
                   hackrf_mode=one(BIT_MODE), qrx=get(BIT_QRX, 4), hrx=one(BIT_HRX), cal=one(BIT_CAL), lna2=one(BIT_LNA2))

    def to_bytes(self) -> bytes:
        return self.to_int().to_bytes(WORD_BYTES, "little")

    @classmethod
    def from_bytes(cls, b: bytes) -> "ControlWord":
        return cls.from_int(int.from_bytes(b[:WORD_BYTES], "little"))


def direction(theta_deg: float, phi_deg: float) -> tuple[float, float]:
    """Direction cosines (u, v) for an angle theta off boresight at azimuth phi (0 = +x, 90 = +y)."""
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return float(np.sin(t) * np.cos(p)), float(np.sin(t) * np.sin(p))


def angles(u: float, v: float) -> tuple[float, float]:
    """Inverse of direction(): (theta_deg, phi_deg)."""
    return float(np.degrees(np.arcsin(min(1.0, np.hypot(u, v))))), float(np.degrees(np.arctan2(v, u)))


def steering(u: float, v: float, f: float = F0) -> np.ndarray:
    """Phase of a plane wave from (u, v) at each element, relative to the array centre."""
    return np.exp(2j * np.pi * f / C0 * (X * u + Y * v))
