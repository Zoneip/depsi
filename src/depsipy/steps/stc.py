"""Spatio-temporal consistency (STC) of PS time series (step8/9).

Port of ``ps_spatio_temporal_consistency.m`` / the official Python
``compute_spatiotemporal_consistency``: for each PS, compare its unwrapped
phase with neighbours in a distance annulus, take the temporal double
differences, and use the minimum std as the STC value.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

__all__ = ["spatiotemporal_consistency"]


def spatiotemporal_consistency(
    phase: np.ndarray,
    az: np.ndarray,
    r: np.ndarray,
    *,
    min_dist: float,
    max_dist: float,
    az_spacing: float = 1.0,
    r_spacing: float = 1.0,
    to_mm: bool = True,
) -> np.ndarray:
    """Return the STC value per PS (``NaN`` when no neighbour is available)."""
    phase = np.atleast_2d(np.asarray(phase, dtype=float))
    az = np.asarray(az, dtype=float).ravel()
    r = np.asarray(r, dtype=float).ravel()
    if phase.shape[0] != az.size or az.size != r.size:
        raise ValueError("phase rows, az and r must have the same size")
    if min_dist < 0 or max_dist < min_dist:
        raise ValueError("require 0 <= min_dist <= max_dist")
    xy = np.column_stack([az * az_spacing, r * r_spacing])
    tree = cKDTree(xy)
    neighbours = tree.query_ball_point(xy, max_dist)
    out = np.full(az.size, np.nan)
    for i, idx in enumerate(neighbours):
        idx = np.asarray(idx, dtype=np.int64)
        idx = idx[(idx != i)]
        if idx.size == 0:
            continue
        dist = np.hypot(xy[idx, 0] - xy[i, 0], xy[idx, 1] - xy[i, 1])
        idx = idx[(dist >= min_dist) & (dist <= max_dist)]
        if idx.size == 0:
            continue
        diff = phase[idx] - phase[i]
        dd = diff[:, :-1] - diff[:, 1:]
        if dd.shape[1] == 0:
            continue
        # neighbours whose time series contain gaps do not vote (``ps_results_*``
        # only holds complete series in the MATLAB chain, but the Python pipeline
        # keeps the removed points as NaN rows)
        complete = np.isfinite(dd).all(axis=1)
        if not complete.any():
            continue
        std_ts = np.std(dd[complete], axis=1)
        if std_ts.size == 0:
            continue
        stc = float(std_ts.min())
        out[i] = stc * 1000.0 if to_mm else stc
    return out
