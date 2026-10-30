"""Orbit interpolation for radar geocoding (port of ``intrp_orbit.m``).

Doris stores the precise orbit as ``[time x y z]`` rows in ECEF.  Geocoding
needs the satellite position, velocity and acceleration at an arbitrary azimuth
time, so the reference implementation fits a degree-2 polynomial per axis in
normalised time, converts the coefficients back to the original time axis, and
differentiates analytically.  The *sampled* state vectors are then spline
interpolated at the requested epochs (MATLAB's ``interp1(..., 'spline')``).

The same module also holds the linear interpolation used by the zero-Doppler
iteration of ``xyz2t.m`` (MATLAB calls ``interp1`` without a method there, i.e.
linear interpolation).
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

__all__ = ["parabolic_state", "interpolate_state", "sample_state_linear"]


def parabolic_state(orbit: np.ndarray) -> np.ndarray:
    """``[t x y z vx vy vz ax ay az]`` sampled at the orbit epochs.

    ``orbit`` is ``(N, 4)``: ``[time, x, y, z]``.  The acceleration of a
    degree-2 polynomial is constant, as in the reference (``repmat(polyder(...))``).
    """
    orbit = np.asarray(orbit, dtype=float)
    if orbit.ndim != 2 or orbit.shape[1] != 4:
        raise ValueError("orbit must be (N, 4) [time, x, y, z]")
    t = orbit[:, 0]
    if t.size < 2:
        raise ValueError("the orbit needs at least two epochs")
    f = float(t.min())
    g = float(t.max() - f)
    if g == 0:
        raise ValueError("the orbit epochs must differ")
    px = (t - f) / g

    pos = [orbit[:, axis] for axis in (1, 2, 3)]
    vel = []
    acc = []
    for axis in range(1, 4):
        c, b, a = np.polyfit(px, orbit[:, axis], 2)  # MATLAB order [c b a]
        coef = np.array(
            [c / g**2, b / g - (2.0 * c * f) / g**2, a - b * f / g + c * (f / g) ** 2]
        )
        vel.append(np.polyval(np.polyder(coef), t))
        acc.append(np.full(t.size, 2.0 * coef[0]))
    # MATLAB order: [t x y z vx vy vz ax ay az]
    return np.column_stack([t, *pos, *vel, *acc])


def interpolate_state(line, orbit: np.ndarray, *, prf: float, ta1: float) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(norm_orbit, norm_orbit_line)`` for the given (1-based) lines.

    ``norm_orbit`` is the sampled state (see :func:`parabolic_state`);
    ``norm_orbit_line`` has one row per input line with
    ``[azimuth_time, x, y, z, vx, vy, vz, ax, ay, az]``.
    """
    norm = parabolic_state(orbit)
    line = np.atleast_1d(np.asarray(line, dtype=float))
    aztime = (line - 1.0) / float(prf) + float(ta1)
    t = norm[:, 0]
    if aztime.min() < t[0] or aztime.max() > t[-1]:
        raise ValueError(
            f"azimuth time(s) {aztime.min():.3f}..{aztime.max():.3f} outside the "
            f"orbit span {t[0]:.3f}..{t[-1]:.3f}"
        )
    pos = CubicSpline(t, norm[:, 1:4])(aztime)
    vel = CubicSpline(t, norm[:, 4:7])(aztime)
    acc = CubicSpline(t, norm[:, 7:10])(aztime)
    out = np.empty((line.size, 10))
    out[:, 0] = aztime
    out[:, 1:4] = np.atleast_2d(pos)
    out[:, 4:7] = np.atleast_2d(vel)
    out[:, 7:10] = np.atleast_2d(acc)
    return norm, out


def sample_state_linear(norm: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Position/velocity/acceleration at ``t`` with linear interpolation.

    Mirrors ``xyz2t.m``, which calls ``interp1`` *without* a method (linear)
    inside the iteration, while using the spline for the final position.
    """
    pos = np.array([np.interp(t, norm[:, 0], norm[:, k]) for k in range(1, 4)])
    vel = np.array([np.interp(t, norm[:, 0], norm[:, k]) for k in range(4, 7)])
    acc = np.array([np.interp(t, norm[:, 0], norm[:, k]) for k in range(7, 10)])
    return pos, vel, acc
