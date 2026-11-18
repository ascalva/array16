"""Array types and physical constants (SI)."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numpy.typing import NDArray

Arr = NDArray[np.float64]
CArr = NDArray[np.complex128]
IntArr = NDArray[np.intp]
BoolArr = NDArray[np.bool_]
Winding = Literal["single", "return"]

E = 1.602176634e-19
ME = 9.1093837015e-31
C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7 * math.pi
KB = 1.380649e-23
Z0 = MU0 * C0
QM_E = -E / ME
