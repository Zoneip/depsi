"""Tests for noise filtering and STC (F3/W4)."""

from __future__ import annotations

import numpy as np

from depsipy.steps.filtering import filter_weights
from depsipy.steps.noise_filter import filter_deformation
from depsipy.steps.stc import spatiotemporal_consistency


def _ref_noise_filter(defo, Btemp, length, kind="gaussian"):
    W = filter_weights(Btemp, length, kind)
    return defo @ W.T


def test_noise_filter_matches_weight_product():
    rng = np.random.default_rng(7)
    Btemp = np.linspace(-1, 1, 12)
    defo = rng.normal(size=(4, 12))
    out = filter_deformation(defo, Btemp, length=1.0)
    np.testing.assert_allclose(out, _ref_noise_filter(defo, Btemp, 1.0), rtol=0, atol=1e-12)
    # smoothing reduces the high-frequency variance
    noisy = np.sin(20 * Btemp)
    assert np.std(filter_deformation(noisy[None, :], Btemp, length=1.0)) < np.std(noisy)


def _ref_stc(phase, az, r, min_dist, max_dist, az_spacing=1.0, r_spacing=1.0):
    xy = np.column_stack([az * az_spacing, r * r_spacing])
    out = np.full(az.size, np.nan)
    for i in range(az.size):
        d = np.hypot(xy[:, 0] - xy[i, 0], xy[:, 1] - xy[i, 1])
        idx = np.nonzero((d >= min_dist) & (d <= max_dist) & (np.arange(az.size) != i))[0]
        if idx.size == 0:
            continue
        diff = phase[idx] - phase[i]
        dd = diff[:, :-1] - diff[:, 1:]
        if dd.shape[1] == 0:
            continue
        out[i] = np.min(np.std(dd, axis=1)) * 1000.0
    return out


def test_stc_matches_bruteforce():
    rng = np.random.default_rng(8)
    n = 20
    az = rng.uniform(0, 20, n)
    r = rng.uniform(0, 20, n)
    phase = rng.normal(size=(n, 5))
    got = spatiotemporal_consistency(phase, az, r, min_dist=1.0, max_dist=6.0)
    ref = _ref_stc(phase, az, r, 1.0, 6.0)
    np.testing.assert_allclose(got, ref, equal_nan=True)


def test_stc_without_neighbours_is_nan():
    out = spatiotemporal_consistency(
        np.zeros((3, 4)), np.array([0, 100, 200]), np.zeros(3),
        min_dist=1.0, max_dist=5.0,
    )
    assert np.isnan(out).all()
