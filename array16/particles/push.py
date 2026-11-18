"""The relativistic Boris step."""

from __future__ import annotations

import math

import numpy as np

from .base import C0, Arr


def cross(a: Arr, b: Arr) -> Arr:
    """Cross product along the last axis."""
    return np.stack(
        [
            a[..., 1] * b[..., 2] - a[..., 2] * b[..., 1],
            a[..., 2] * b[..., 0] - a[..., 0] * b[..., 2],
            a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0],
        ],
        axis=-1,
    )


def push(
    x: Arr, u: Arr, qm: Arr, dt: float, ef: Arr, bf: Arr, drag: float = 0.0, acc: Arr | None = None
) -> tuple[Arr, Arr]:
    """One relativistic Boris step for N particles: positions x (N, 3), u = gamma v (N, 3), charge/mass qm (N,).

    The magnetic rotation uses the tangent of the exact half angle, so the gyration frequency has no time-step
    error. The price is in crossed fields: the E x B drift comes out low by (theta/2) / tan(theta/2), theta being
    the gyro-angle per step (0.2 % at 40 steps per turn). drag is a Stokes rate (1/s); acc is a uniform
    acceleration (3,) in m/s^2.
    """
    k = (qm * dt / 2)[:, None]
    um = u + k * ef
    if acc is not None:
        um = um + acc * (dt / 2)
    gam = np.sqrt(1 + np.sum(um * um, 1) / C0**2)
    bn = np.sqrt(np.sum(bf * bf, 1))
    tv = bf * (np.tan(qm * dt / 2 * bn / gam) / np.where(bn > 0, bn, 1.0))[:, None]
    up = um + cross(um, tv)
    u = um + cross(up, tv * (2 / (1 + np.sum(tv * tv, 1)))[:, None]) + k * ef
    if acc is not None:
        u = u + acc * (dt / 2)
    if drag:
        u = u * math.exp(-drag * dt)
    gam = np.sqrt(1 + np.sum(u * u, 1) / C0**2)
    return x + u / gam[:, None] * dt, u
