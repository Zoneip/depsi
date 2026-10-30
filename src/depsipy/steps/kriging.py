"""Variogram fitting and ordinary Kriging (step4).

The covariance/variogram conventions follow ``ps_kriging.m``/``ps_variogram.m``:
``gamma(h) = c1 * shape(h/a) + c2`` and ``C(h) = c1 * corr(h/a) + c2 * I``.
Only the numerically stable ordinary-Kriging path is implemented here.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

__all__ = ["variogram_cloud", "fit_variogram", "kriging_predict"]


def _shape(h: np.ndarray, a: float, model: str) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    if a <= 0:
        raise ValueError("variogram range must be positive")
    if model == "exponential":
        return 1.0 - np.exp(-h / a)
    if model == "gaussian":
        return 1.0 - np.exp(-(h / a) ** 2)
    if model == "spherical":
        u = np.minimum(h / a, 1.0)
        return 1.5 * u - 0.5 * u ** 3
    raise ValueError(f"unknown variogram model {model!r}")


def _corr(h: np.ndarray, a: float, model: str) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    if model == "exponential":
        return np.exp(-h / a)
    if model == "gaussian":
        return np.exp(-(h / a) ** 2)
    if model == "spherical":
        u = np.minimum(h / a, 1.0)
        return 1.0 - 1.5 * u + 0.5 * u ** 3
    raise ValueError(f"unknown variogram model {model!r}")


def variogram_cloud(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    cutoff: float = 10000.0,
    nlags: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    """Pairwise distance / semivariance bins (an omnidirectional cloud)."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    if x.size != y.size or x.size != z.size:
        raise ValueError("x, y and z must have the same size")
    if x.size < 2:
        raise ValueError("at least two points are needed")
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    dz = z[:, None] - z[None, :]
    dist = np.sqrt(dx * dx + dy * dy)
    semi = 0.5 * dz * dz
    upper = np.triu(np.ones((x.size, x.size), dtype=bool), 1)
    keep = upper & (dist <= cutoff)
    d, s = dist[keep], semi[keep]
    if d.size == 0:
        raise ValueError("no pairs within cutoff")
    bins = np.linspace(0.0, cutoff, int(nlags) + 1)
    which = np.digitize(d, bins, right=True) - 1
    which = np.clip(which, 0, int(nlags) - 1)
    lags = np.array([d[which == i].mean() if np.any(which == i) else np.nan for i in range(int(nlags))])
    vals = np.array([s[which == i].mean() if np.any(which == i) else np.nan for i in range(int(nlags))])
    valid = np.isfinite(lags) & np.isfinite(vals)
    return lags[valid], vals[valid]


def fit_variogram(
    lags: np.ndarray,
    semivariance: np.ndarray,
    model: str = "gaussian",
) -> tuple[float, float, float]:
    """Fit ``gamma(h)=c1*shape(h/a)+c2``; returns ``(c1, a, c2)``."""
    lags = np.asarray(lags, dtype=float)
    semi = np.asarray(semivariance, dtype=float)
    if lags.size < 3 or lags.size != semi.size:
        raise ValueError("need at least three lag/semivariance pairs")
    if model not in ("exponential", "gaussian", "spherical"):
        raise ValueError(f"unknown variogram model {model!r}")

    def f(h, c1, a, c2):
        return c1 * _shape(h, a, model) + c2

    c1_0 = max(float(np.nanmax(semi)), 1e-12)
    a_0 = max(float(np.nanmax(lags)) / 3.0, 1e-9)
    try:
        popt, _ = curve_fit(
            f, lags, semi,
            p0=[c1_0, a_0, 0.0],
            bounds=([0.0, 1e-12, 0.0], [np.inf, np.inf, np.inf]),
            maxfev=10000,
        )
        c1, a, c2 = (float(v) for v in popt)
        if not np.all(np.isfinite([c1, a, c2])) or a <= 0:
            raise RuntimeError("invalid fit")
        return c1, a, c2
    except Exception:
        return c1_0, a_0, 0.0


def kriging_predict(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    x0: np.ndarray,
    y0: np.ndarray,
    *,
    model: str = "gaussian",
    c1: float,
    a: float,
    c2: float = 0.0,
    nmax: int = 16,
    max_dist: float = np.inf,
    method: str = "ordinary",
) -> np.ndarray:
    """Predict ``z`` at ``(x0,y0)`` with local ordinary Kriging."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    x0 = np.asarray(x0, dtype=float).ravel()
    y0 = np.asarray(y0, dtype=float).ravel()
    if x.size != y.size or x.size != z.size:
        raise ValueError("x, y and z must have the same size")
    if x0.size != y0.size:
        raise ValueError("x0 and y0 must have the same size")
    out = np.full(x0.size, np.nan)
    n = min(int(nmax), x.size)
    if n < 1 or x.size == 0:
        return out

    for k in range(x0.size):
        d = np.hypot(x - x0[k], y - y0[k])
        idx = np.argsort(d)[:n]
        idx = idx[d[idx] <= max_dist]
        if idx.size < 3:
            continue
        xs, ys, zs = x[idx], y[idx], z[idx]
        h = np.hypot(xs[:, None] - xs[None, :], ys[:, None] - ys[None, :])
        kxx = c1 * _corr(h, a, model)
        kxx[np.diag_indices_from(kxx)] += c2
        h0 = np.hypot(xs - x0[k], ys - y0[k])
        kx0 = c1 * _corr(h0, a, model)
        mean_z = zs.mean()
        rhs = zs - mean_z

        if method == "simple":
            weights = np.linalg.solve(kxx, kx0)
            out[k] = float(weights @ rhs + mean_z)
            continue
        if method != "ordinary":
            raise ValueError(f"unknown kriging method {method!r}")
        m = idx.size
        system = np.empty((m + 1, m + 1))
        system[:m, :m] = kxx
        system[:m, m] = 1.0
        system[m, :m] = 1.0
        system[m, m] = 0.0
        target = np.empty(m + 1)
        target[:m] = kx0
        target[m] = 1.0
        try:
            weights = np.linalg.solve(system, target)[:m]
        except np.linalg.LinAlgError:
            continue
        out[k] = float(weights @ rhs + mean_z)
    return out
