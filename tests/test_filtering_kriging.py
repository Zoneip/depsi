"""Tests for temporal filtering and Kriging (F2/W3)."""

from __future__ import annotations

import numpy as np

from depsipy.steps.filtering import filter_weights, separate_atmosphere
from depsipy.steps.kriging import fit_variogram, kriging_predict, variogram_cloud


def _ref_weights(Btemp, length, kind):
    span = float(Btemp.max() - Btemp.min())
    n = int(np.ceil(span * 1000)) + 1
    if kind == "block":
        kernel = np.zeros(n)
        kernel[: int(round(0.5 * length * 1000)) + 1] = 1.0
    elif kind == "triangle":
        half = int(round(0.5 * length * 1000))
        kernel = np.zeros(n)
        kernel[: half + 1] = np.linspace(0, 1, half + 1)[::-1]
    elif kind == "gaussian":
        sigma = 0.5 * (length / 3.0) * 1000.0
        kernel = np.exp(-0.5 * (np.arange(n) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)
    W = np.zeros((Btemp.size, Btemp.size))
    for v in range(Btemp.size):
        idx = np.floor(np.abs(Btemp - Btemp[v]) * 1000).astype(int)
        idx = np.clip(idx, 0, n - 1)
        W[v] = kernel[idx]
        W[v] /= W[v].sum()
    return W


def test_filter_weights_match_loop():
    Btemp = np.array([-0.8, -0.3, 0.0, 0.4, 0.9])
    for kind in ("block", "triangle", "gaussian"):
        got = filter_weights(Btemp, 1.5, kind)
        np.testing.assert_allclose(got, _ref_weights(Btemp, 1.5, kind), rtol=0, atol=1e-12)


def test_gaussian_filter_separates_linear_and_high_frequency():
    Btemp = np.linspace(-1.0, 1.0, 25)
    trend = 2.0 * Btemp
    high = 0.3 * np.sin(12 * np.pi * Btemp)
    phase = trend + high
    non_linear, atmosphere = separate_atmosphere(phase[None, :], Btemp, length=0.8)
    # low-pass keeps the trend and removes most of the fast oscillation;
    # the first/last samples are boundary-affected, so check the interior.
    np.testing.assert_allclose(non_linear[0, 2:-2], trend[2:-2], atol=0.05)
    assert np.std(atmosphere[0]) < 0.2


def test_block_filter_and_master_atmo():
    Btemp = np.array([0.0, 0.2, 0.4])
    phase = np.array([[1.0, 2.0, 3.0]])
    non_linear, atmo = separate_atmosphere(phase, Btemp, length=1.0, kind="block", master_atmo=np.array([0.5]))
    np.testing.assert_allclose(non_linear, [[2.0, 2.0, 2.0]], atol=1e-12)
    np.testing.assert_allclose(atmo, [[-0.5, 0.5, 1.5]], atol=1e-12)


def test_variogram_cloud():
    x = np.array([0.0, 1.0, 2.0])
    y = np.zeros(3)
    z = np.array([0.0, 1.0, 2.0])
    lags, semis = variogram_cloud(x, y, z, cutoff=2.0, nlags=2)
    np.testing.assert_allclose(lags, [1.0, 2.0])
    np.testing.assert_allclose(semis, [0.5, 2.0])


def test_kriging_constant_field():
    x = np.array([0.0, 1.0, 0.0, 1.0, 0.5])
    y = np.array([0.0, 0.0, 1.0, 1.0, 0.5])
    z = np.full(5, 5.0)
    out = kriging_predict(x, y, z, np.array([0.2, 0.8]), np.array([0.3, 0.7]),
                          model="gaussian", c1=1.0, a=1.0, c2=0.0, nmax=5)
    np.testing.assert_allclose(out, 5.0, atol=1e-8)


def test_kriging_smooth_field_is_close():
    rng = np.random.default_rng(2)
    x = rng.uniform(0, 10, 40)
    y = rng.uniform(0, 10, 40)
    z = np.sin(x / 3.0) + np.cos(y / 4.0)
    x0 = rng.uniform(1, 9, 10)
    y0 = rng.uniform(1, 9, 10)
    lags, semis = variogram_cloud(x, y, z, cutoff=20.0, nlags=25)
    c1, a, c2 = fit_variogram(lags, semis, "gaussian")
    pred = kriging_predict(x, y, z, x0, y0, model="gaussian", c1=c1, a=a, c2=c2, nmax=12)
    truth = np.sin(x0 / 3.0) + np.cos(y0 / 4.0)
    assert np.isfinite(pred).all()
    assert np.sqrt(np.mean((pred - truth) ** 2)) < 0.5
