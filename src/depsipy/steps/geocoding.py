"""Geodetic conversions and radar geocoding for output (step7).

``latlonh_to_xyz``/``xyz_to_latlonh`` are the WGS84 ellipsoid <-> ECEF
conversions.  On top of them, the radar geometry of a Doris result file is
implemented here:

* :func:`lph2xyz` -- ``(line, pixel, height)`` -> ECEF, solving zero-Doppler,
  range and ellipsoid equations (port of ``lph2xyz.m``);
* :func:`xyz2t` / :func:`xyz2lp` -- the inverse direction (port of
  ``xyz2t.m`` / ``xyz2lp.m``);
* :func:`radar_to_geodetic` -- the convenient ``(line, pixel, height)`` ->
  ``(lat, lon, h)`` combination used by ``ps_output.m``.

The metadata container and the precise orbit come from
:func:`depsipy.io.doris.read_doris_metadata`; the orbit interpolation lives in
:mod:`depsipy.core.orbit`.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from ..core.orbit import interpolate_state, sample_state_linear

__all__ = [
    "latlonh_to_xyz",
    "xyz_to_latlonh",
    "orbit_state",
    "lph2xyz",
    "xyz2t",
    "xyz2lp",
    "radar_to_geodetic",
]

# WGS84
_A = 6378137.0
_F = 1.0 / 298.257223563
_E2 = _F * (2.0 - _F)
_B = _A * (1.0 - _F)

_SOL = 299792458.0  # speed of light [m/s]


def latlonh_to_xyz(lat: np.ndarray, lon: np.ndarray, height: np.ndarray) -> np.ndarray:
    """WGS84 geodetic (rad, rad, m) -> ECEF (m), shape ``(..., 3)``."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    h = np.asarray(height, dtype=float)
    lat, lon, h = np.broadcast_arrays(lat, lon, h)
    n = _A / np.sqrt(1.0 - _E2 * np.sin(lat) ** 2)
    x = (n + h) * np.cos(lat) * np.cos(lon)
    y = (n + h) * np.cos(lat) * np.sin(lon)
    z = (n * (1.0 - _E2) + h) * np.sin(lat)
    return np.stack([x, y, z], axis=-1)


def xyz_to_latlonh(xyz: np.ndarray, *, iterations: int = 12) -> np.ndarray:
    """ECEF (m) -> WGS84 geodetic (rad, rad, m), shape ``(..., 3)``."""
    xyz = np.asarray(xyz, dtype=float)
    if xyz.shape[-1] != 3:
        raise ValueError("last dimension must be 3")
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1.0 - _E2))
    for _ in range(int(iterations)):
        n = _A / np.sqrt(1.0 - _E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1.0 - _E2 * n / (n + h)))
    n = _A / np.sqrt(1.0 - _E2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - n
    return np.stack([lat, lon, h], axis=-1)


# ---------------------------------------------------------------------------
# radar geometry (Doris result file)
# ---------------------------------------------------------------------------


def orbit_state(meta, line) -> tuple[np.ndarray, np.ndarray]:
    """``intrp_orbit.m`` for a Doris metadata container."""
    return interpolate_state(line, meta.orbit, prf=meta.prf, ta1=meta.ta1)


def lph2xyz(
    line,
    pixel,
    height,
    meta,
    orbit_line: np.ndarray,
    *,
    max_iter: int = 100,
    crit: float = 1e-8,
) -> np.ndarray:
    """``(line, pixel, height)`` -> ECEF for a Doris product (``lph2xyz.m``).

    Solves the zero-Doppler, range and ellipsoid equations simultaneously with a
    Newton iteration, initialised at the scene centre (as in the reference).
    ``orbit_line`` is the output of :func:`orbit_state`.
    """
    line = np.atleast_1d(np.asarray(line, dtype=float))
    pixel = np.atleast_1d(np.asarray(pixel, dtype=float))
    height = np.atleast_1d(np.asarray(height, dtype=float))
    if line.size != pixel.size or line.size != height.size:
        raise ValueError("line, pixel and height must have the same size")
    if orbit_line.shape[0] != line.size:
        raise ValueError("orbit_line must have one row per point")

    # scene-centre ECEF, used as the initial guess
    h0 = float(np.mean(height))
    phi = np.deg2rad(meta.center_lat)
    lam = np.deg2rad(meta.center_lon)
    n_center = _A / np.sqrt(1.0 - _E2 * np.sin(phi) ** 2)
    scene = np.array(
        [
            (n_center + h0) * np.cos(phi) * np.cos(lam),
            (n_center + h0) * np.cos(phi) * np.sin(lam),
            (n_center * (1.0 - _E2) + h0) * np.sin(phi),
        ]
    )

    out = np.empty((line.size, 3))
    for k in range(line.size):
        pos = scene.copy()
        sat = orbit_line[k, 1:4]
        vel = orbit_line[k, 4:7]
        ratime = meta.tr1 + (pixel[k] - 1.0) / meta.rsr
        a_h = _A + height[k]
        b_h = _B + height[k]
        for _ in range(int(max_iter)):
            d = pos - sat
            rhs = np.array(
                [
                    -(vel @ d),
                    -(d @ d - (_SOL * ratime) ** 2),
                    -((pos[0] ** 2 + pos[1] ** 2) / a_h**2 + (pos[2] / b_h) ** 2 - 1.0),
                ]
            )
            jac = np.array(
                [
                    [vel[0], vel[1], vel[2]],
                    [2.0 * d[0], 2.0 * d[1], 2.0 * d[2]],
                    [2.0 * pos[0] / a_h**2, 2.0 * pos[1] / a_h**2, 2.0 * pos[2] / b_h**2],
                ]
            )
            step = np.linalg.solve(jac, rhs)
            pos = pos + step
            if np.all(np.abs(step) < crit):
                break
        else:  # pragma: no cover - only for badly conditioned geometry
            raise RuntimeError(f"lph2xyz did not converge for point {k} (line {line[k]})")
        out[k] = pos
    return out


def xyz2t(xyz: np.ndarray, meta, *, max_iter: int = 100, crit: float = 1e-8) -> tuple[np.ndarray, np.ndarray]:
    """ECEF -> (azimuth time, range time) with zero-Doppler (``xyz2t.m``)."""
    xyz = np.atleast_2d(np.asarray(xyz, dtype=float))
    if xyz.shape[1] != 3:
        raise ValueError("xyz must be (N, 3)")
    norm, _ = orbit_state(meta, meta.mid_line)
    ta1 = float(meta.ta1)
    prf = float(meta.prf)
    start = (float(meta.mid_line) - 1.0) / prf + ta1

    t_azi = np.empty(xyz.shape[0])
    t_ran = np.empty(xyz.shape[0])
    for k in range(xyz.shape[0]):
        pos = xyz[k]
        tazi = start
        for _ in range(int(max_iter)):
            sat, vel, acc = sample_state_linear(norm, tazi)
            delta = pos - sat
            step = -(vel @ delta) / (acc @ delta - vel @ vel)
            tazi += step
            if abs(step) < crit:
                break
        t_azi[k] = tazi
        # the reference uses the spline for the final satellite position
        sat_final = np.array([float(CubicSpline(norm[:, 0], norm[:, j])(tazi)) for j in (1, 2, 3)])
        t_ran[k] = np.linalg.norm(pos - sat_final) / _SOL
    return t_azi, t_ran


def xyz2lp(xyz: np.ndarray, meta, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """ECEF -> (line, pixel) for a Doris product (``xyz2lp.m``)."""
    t_azi, t_ran = xyz2t(xyz, meta, **kwargs)
    line = 1.0 + meta.prf * (t_azi - meta.ta1)
    pixel = 1.0 + meta.rsr * (t_ran - meta.tr1)
    return line, pixel


def radar_to_geodetic(line, pixel, height, meta) -> np.ndarray:
    """``(line, pixel, height)`` -> ``(lat, lon, h)`` in degrees (``ps_output.m``)."""
    _, orbit_line = orbit_state(meta, line)
    xyz = lph2xyz(line, pixel, height, meta, orbit_line)
    latlonh = xyz_to_latlonh(xyz)
    latlonh[:, :2] = np.rad2deg(latlonh[:, :2])
    return latlonh
