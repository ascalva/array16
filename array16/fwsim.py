"""Reference model of the beam-controller firmware (RP2040).

It consumes the same frames as the real firmware will and answers the one
question the simulator needs: which control word is latched at time t.
The real firmware should behave exactly like this model; that is the spec.
tests/test_firmware_core.py runs the C core and this model on the same frames, and tests/test_firmware_emulated.py the whole
firmware in the rp2040js emulator.

Sequences: RUN latches a blank marker word (all zeros) at once, then the
loaded words one slot each, and repeats. With repeat = 0 that goes on until
STOP, SET_WORD or LOAD_SEQ. With repeat = n it ends after n periods, and the
controller then behaves as after STOP: the static word (the last SET_WORD)
is latched again and STATUS reports that nothing runs.

Frames are protocol version 2 (p.FW_VERSION): LOAD_SEQ is slot_us u16, repeat u16,
count u16, then count words of 7 bytes, and is refused with ERR_ARG unless the count
is 1 to p.MAX_SEQ = 256, the slot at least 10 us and the payload exactly 6 + 7 * count
bytes; STATUS answers with the 10 bytes of p.Status, its count u16.
"""
from __future__ import annotations

import struct

import numpy as np

from . import protocol as p
from .board import ControlWord


class ControllerModel:
    def __init__(self):
        self.word = ControlWord()          # static word, used when no sequence runs
        self.seq: p.Sequence | None = None
        self.running = False
        self.t_start = 0.0                 # time the running sequence started, s
        self.clock_src = 0                 # 0 = on-board 40 MHz, 1 = HackRF 10 MHz
        self.ref_ok = True                 # the selected reference is present and times the slots
        self.usb_ma = 3000                 # current the USB-C source advertises, mA
        self.now = 0.0                     # set by the backend before each frame

    # ---- frame handling ------------------------------------------------------
    def handle(self, frame: bytes) -> bytes:
        if self.running and self._finished(self.now):
            self.running = False
        parsed = p.parse_request(frame)
        if parsed is None:
            return p.reply(0, p.ERR_CRC)
        cmd, payload = parsed
        try:
            if cmd == p.PING:
                return p.reply(cmd, p.OK, bytes([p.FW_VERSION]))
            if cmd == p.SET_WORD:
                if len(payload) != 7:
                    raise ValueError
                self.word, self.running = ControlWord.from_bytes(payload), False
                return p.reply(cmd)
            if cmd == p.LOAD_SEQ:
                self.seq, self.running = p.Sequence.from_payload(payload), False
                return p.reply(cmd)
            if cmd == p.RUN:
                if self.seq is None:
                    raise ValueError
                self.running, self.t_start = True, self.now
                return p.reply(cmd)
            if cmd == p.STOP:
                self.running = False
                return p.reply(cmd)
            if cmd == p.CLOCK_SRC:
                if len(payload) != 1 or payload[0] > 1:
                    raise ValueError
                self.clock_src = payload[0]
                return p.reply(cmd)
            if cmd == p.STATUS:
                return p.reply(cmd, p.OK, p.Status(p.FW_VERSION, self.running, len(self.seq.words) if self.seq else 0,
                                                   self.seq.slot_us if self.seq else 0, self.clock_src, self.ref_ok,
                                                   self.usb_ma).payload())
        except (ValueError, struct.error):
            return p.reply(cmd, p.ERR_ARG)
        return p.reply(cmd, p.ERR_CMD)

    # ---- timing --------------------------------------------------------------
    def _finished(self, t: float) -> bool:
        if not (self.running and self.seq is not None and self.seq.repeat):
            return False
        period = (len(self.seq.words) + 1) * self.seq.slot_us * 1e-6
        return t - self.t_start >= self.seq.repeat * period - 1e-12

    def words(self) -> list[ControlWord]:
        """Word table for the current state.

        While a sequence runs: the blank marker word, the sequence, then the
        static word that returns when a finite sequence ends. Otherwise just the static word.
        """
        if self.running and self.seq is not None:
            return [ControlWord()] + list(self.seq.words) + [self.word]
        return [self.word]

    def word_index(self, t: np.ndarray) -> np.ndarray:
        """Index into words() of the word latched at each time in t (seconds)."""
        if not (self.running and self.seq is not None):
            return np.zeros(len(t), dtype=int)
        slot = self.seq.slot_us * 1e-6
        k = np.floor((t - self.t_start) / slot + 1e-9).astype(int)
        period = len(self.seq.words) + 1
        idx = k % period
        after = period                                   # index of the static word in words()
        if self.seq.repeat:
            idx = np.where(k >= self.seq.repeat * period, after, idx)
        return np.where(k < 0, after, idx)
