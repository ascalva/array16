"""What api.Array needs from a backend. sim.SimBackend implements it; the hardware backend will too.

On hardware a session will open the port, then call api.Array.check_version: the controller must speak
protocol.FW_VERSION.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


@dataclass
class Capture:
    samples: np.ndarray           # (channels, n) complex64, full scale = 1
    fs: float
    scale: float                  # sqrt(W) at the radio input per unit sample


class Backend(Protocol):
    kind: str                     # "quad" (four channels on two radio boards) or "hackrf" (one channel)
    fs: float                     # sample rate, Hz

    def controller(self, frame: bytes) -> bytes:
        """Send one request frame to the beam controller and return its reply frame."""

    def tune(self, f: float) -> None:
        """Retune every radio channel to f (Hz). Receive and transmit oscillator phases change."""

    def capture(self, n: int) -> Capture:
        """n samples per channel, taken while the controller does whatever it was last told."""

    def cal_signal(self, on: bool) -> np.ndarray:
        """Four-channel only: start or stop the probe code on the transmit channel cabled to CAL IN.

        Returns one period of the code at the sample rate, the replica to correlate against.
        """

    def set_tx(self, weights: np.ndarray) -> None:
        """Four-channel only: complex baseband weight (|w| <= 1) on each of the four transmit channels."""
