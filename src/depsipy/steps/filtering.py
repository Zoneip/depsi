"""Temporal filtering for atmosphere / non-linear deformation separation (step4).

Ports the ``ps_filtering.m`` filter construction and applies it to a
``(Nps, Nifgs)`` residual matrix.  The MATLAB code loops over epochs; here the
weight matrix is built once and the filtering is one matrix product.
"""

from __future__ import annotations

import numpy as np

__all__ = ["filter_weights", "separate_atmosphere"]


def _filter_kernel(time_span: float, length: float, kind: str) -> np.ndarray:
    if time_span <= 0:
        raise ValueError("time_span must be positive")
    if length <= 0:
        raise ValueError("filter length must be positive")
    n = int(np.ceil(time_span * 1000)) + 1
    kind = str(kind).lower()
    if kind == "block":
        kernel = np.zeros(n)
        kernel[: int(round(0.5 * length * 1000)) + 1] = 1.0
        return kernel
    if kind == "triangle":
        half = int(round(0.5 * length * 1000))
        ramp = np.linspace(0.0, 1.0, half + 1)
        kernel = np.zeros(n)
        kernel[: half + 1] = ramp[::-1]
        return kernel
    if kind == "gaussian":
        sigma = 0.5 * (length / 3.0) * 1000.0
        if sigma <= 0:
            raise ValueError("gaussian filter length is too small")
        h = np.arange(n)
        return np.exp(-0.5 * (h / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)
    raise ValueError(f"unknown filter type {kind!r}")


def filter_weights(
    Btemp: np.ndarray,
    length: float,
    kind: str = "gaussian",
    sampling: int = 1000,
) -> np.ndarray:
    """Return the normalized ``(Nifgs, Nifgs)`` filter weight matrix."""
    Btemp = np.asarray(Btemp, dtype=float)
    if Btemp.ndim != 1 or Btemp.size == 0:
        raise ValueError("Btemp must be a non-empty 1-D array")
    span = float(Btemp.max() - Btemp.min())
    kernel = _filter_kernel(span, float(length), kind)
    diff = np.abs(Btemp[:, None] - Btemp[None, :]) * sampling
    idx = np.clip(np.floor(diff).astype(np.int64), 0, kernel.size - 1)
    weights = kernel[idx]
    total = weights.sum(axis=1, keepdims=True)
    if np.any(total == 0):
        raise ValueError("filter weights sum to zero")
    return weights / total


def separate_atmosphere(
    phase_residuals: np.ndarray,
    Btemp: np.ndarray,
    *,
    length: float = 3.0,
    kind: str = "gaussian",
    master_atmo: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Split residual phase into non-linear deformation and atmosphere.

    Parameters
    ----------
    phase_residuals : ``(Nps, Nifgs)`` residual phase [rad].
    Btemp : ``(Nifgs,)`` temporal baselines [year].
    master_atmo : optional ``(Nps,)`` master-atmosphere phase to add back.

    Returns
    -------
    ``(non_linear, atmosphere)`` with the same shape as ``phase_residuals``.
    """
    phase = np.atleast_2d(np.asarray(phase_residuals, dtype=float))
    Btemp = np.asarray(Btemp, dtype=float)
    if phase.shape[1] != Btemp.size:
        raise ValueError("phase_residuals second dimension must match Btemp")
    weights = filter_weights(Btemp, length, kind)
    non_linear = phase @ weights.T
    atmosphere = phase - non_linear
    if master_atmo is not None:
        atmosphere = atmosphere + np.asarray(master_atmo, dtype=float)[:, None]
    return non_linear, atmosphere
