"""Temporal noise filtering (step8) -- Gaussian smoothing of PS deformation.

The MATLAB ``ps_noise_filter.m`` builds the same Gaussian kernel as
``ps_filtering.m`` and applies it to each PS deformation time series.  Here the
kernel is shared with :mod:`depsipy.steps.filtering` and the filtering is one
matrix product.
"""

from __future__ import annotations

import numpy as np

from .filtering import filter_weights

__all__ = ["filter_deformation"]


def filter_deformation(
    deformation: np.ndarray,
    Btemp: np.ndarray,
    *,
    length: float = 10.0 / 12.0,
    kind: str = "gaussian",
) -> np.ndarray:
    """Smooth ``(Nps, Nifgs)`` deformation time series along time."""
    defo = np.atleast_2d(np.asarray(deformation, dtype=float))
    Btemp = np.asarray(Btemp, dtype=float)
    if defo.shape[1] != Btemp.size:
        raise ValueError("deformation second dimension must match Btemp")
    weights = filter_weights(Btemp, length, kind)
    return defo @ weights.T
