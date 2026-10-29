"""Tests for the LAMBDA integer least-squares / bootstrapping core.

Ground truth is *brute-force enumeration* on the integer lattice: for small
dimension we can enumerate every integer vector in a box and compare the squared
norm in the metric of the ambiguity covariance.  That is exactly the quantity
the branch-and-bound search minimises.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from depsipy.core import _kernels as kernels
from depsipy.core._reference import chistart, decorrel, ldldecom, lsearch
from depsipy.core.estimation import build_geometry, evaluate
from depsipy.core.lambda_ils import LambdaModel, resolve_ambiguities
from depsipy.core.models import build_functional_model, model_spec, npar


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def random_spd(n, rng, scale=1.0):
    A = rng.normal(size=(n, n))
    return (A @ A.T + np.eye(n) * 0.5) * scale


def realistic_geometry(n=8, seed=0, model=2, topo_std=30.0):
    """A geometry with perpendicular baselines that actually vary.

    A constant height-to-phase factor makes topography indistinguishable from the
    master atmosphere term, i.e. the functional model becomes rank deficient.
    """
    rng = np.random.default_rng(seed)
    Btemp = np.sort(rng.uniform(-1.2, 1.2, n))
    bperp = np.linspace(-400, 500, n)
    h2ph = -2 * bperp / (0.056 * 800e3 * np.sin(np.deg2rad(23)))
    std = np.array([topo_std, 0.005, 1, 0.02, 0.01, 0.01, 0.005])
    return build_geometry(
        model, n, h2ph, Btemp, np.zeros(n), std, None,
        m2ph=-4 * np.pi / 0.056, sat_vel=7550.0,
    )


# ---------------------------------------------------------------------------
# LtDL decomposition / decorrelation
# ---------------------------------------------------------------------------


def test_ldldecom_reconstructs_matrix():
    rng = np.random.default_rng(0)
    Q = random_spd(6, rng)
    L, D = ldldecom(Q)
    # DePSI/LAMBDA eliminates from the last variable, so Q = L.T diag(D) L
    np.testing.assert_allclose(L.T @ np.diag(D) @ L, Q, atol=1e-10)
    np.testing.assert_allclose(np.diag(L), 1.0)
    np.testing.assert_allclose(np.triu(L, 1), 0.0)


def test_ldldecom_rejects_non_pd():
    with pytest.raises(ValueError):
        ldldecom(np.array([[1.0, 2.0], [2.0, 1.0]]))


def test_decorrel_is_unimodular_and_reduces_correlation():
    rng = np.random.default_rng(1)
    Q = random_spd(5, rng) * 100.0
    Z, L, D = decorrel(Q)
    assert abs(round(np.linalg.det(Z))) == 1
    Zi = np.linalg.inv(Z.astype(float))
    np.testing.assert_allclose(Zi, np.round(Zi), atol=1e-8)
    dec = Z.T @ Q @ Z
    np.testing.assert_allclose(L.T @ np.diag(D) @ L, dec, atol=1e-8)
    # the decorrelated problem is better conditioned than the original
    assert D.max() / D.min() < np.linalg.cond(Q)


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_lsearch_finds_the_global_minimum(seed):
    rng = np.random.default_rng(seed)
    n = 4
    Q = random_spd(n, rng)
    Z, L, D = decorrel(Q)
    Qd = Z.T @ Q @ Z
    Qi = np.linalg.inv(Qd)
    invL, invD = np.linalg.inv(L), 1.0 / D

    afloat = rng.normal(size=n) * 3.0
    best_val = np.inf
    best = None
    for z in itertools.product(range(-8, 9), repeat=n):
        z = np.asarray(z, dtype=float)
        d = afloat - z
        val = d @ Qi @ d
        if val < best_val:
            best_val, best = val, z

    got, sqnorm, ierr, brk = lsearch(afloat, invL, invD, best_val + 0.5)
    assert not brk
    np.testing.assert_array_equal(got, best)
    np.testing.assert_allclose(sqnorm, best_val, rtol=1e-8)


def test_lsearch_returns_nan_when_the_guard_trips():
    rng = np.random.default_rng(2)
    Q = random_spd(4, rng)
    _, L, D = decorrel(Q)
    got, sq, ierr, brk = lsearch(
        rng.normal(size=4) * 2, np.linalg.inv(L), 1.0 / D, 1e4, max_backtrack=5
    )
    assert brk
    assert np.all(np.isnan(got))


def test_chistart_bootstrap_point_is_inside_its_ellipsoid():
    rng = np.random.default_rng(7)
    Q = random_spd(5, rng)
    Z, L, D = decorrel(Q)
    a = rng.normal(size=5)
    chi2, afixed = chistart(L, D, a, 1)
    assert afixed.shape == (5,)
    np.testing.assert_allclose(afixed, np.round(afixed))
    d = a - afixed
    inv = np.linalg.inv(L.T @ np.diag(D) @ L)
    assert d @ inv @ d <= chi2 + 1e-9


@pytest.mark.parametrize("trial", range(12))
def test_numba_and_reference_agree(trial):
    """The compiled kernel must reproduce the reference arithmetic exactly."""
    if not kernels.HAVE_NUMBA:
        pytest.skip("numba not installed")
    rng = np.random.default_rng(trial)
    n = int(rng.integers(3, 9))
    Q = random_spd(n, rng)
    Q = Q @ np.diag(10.0 ** rng.uniform(-2, 2, n)) @ Q.T  # nasty conditioning
    Z, L, D = decorrel(Q)
    invLtDL = np.linalg.inv(L.T @ np.diag(D) @ L)
    a = rng.normal(size=n) * 2

    chi_ref, af_ref = chistart(L, D, a.copy(), 1)
    chi_nb, af_nb = kernels.nb_chistart(L, D, invLtDL, a.copy())
    np.testing.assert_allclose(chi_nb, chi_ref, rtol=1e-12)
    np.testing.assert_array_equal(af_nb, af_ref)

    invL, invD = np.linalg.inv(L), 1.0 / D
    ref = lsearch(a, invL, invD, chi_ref)
    nb = kernels.nb_lsearch(a, invL, invD, chi_ref)
    np.testing.assert_array_equal(nb[0], ref[0])
    np.testing.assert_allclose(nb[1], ref[1], rtol=1e-10)
    assert nb[2] == ref[2]
    assert nb[3] == ref[3]


# ---------------------------------------------------------------------------
# end-to-end estimation on a network geometry
# ---------------------------------------------------------------------------


def test_geometry_recovers_known_parameters():
    """Noise-free arcs must be unwrapped exactly onto the true model."""
    geom = realistic_geometry()
    n = geom.Qahat.shape[0]
    param = np.array([0.002, 0.0, 0.001])  # topo [m], master, linear [m/y]
    phase = geom.B1 @ param
    rng = np.random.default_rng(1)
    amb = rng.integers(-2, 3, (50, n)).astype(float)
    dphase = np.mod(phase[None, :] + 2 * np.pi * amb + np.pi, 2 * np.pi) - np.pi

    res = evaluate(geom, dphase)
    residual = (dphase + 2 * np.pi * res.acheck) - res.param @ geom.B1.T
    np.testing.assert_allclose(residual, 0.0, atol=1e-8)
    np.testing.assert_allclose(res.param, np.tile(param, (50, 1)), atol=1e-6)
    np.testing.assert_allclose(res.ens_coh, 1.0, atol=1e-9)


def test_geometry_recovers_ambiguities_with_informative_prior():
    """A tight height prior makes the integer ambiguity observable; the ILS
    solution must then be exact (this is the regime the estimator relies on)."""
    geom = realistic_geometry(topo_std=0.5)
    n = geom.Qahat.shape[0]
    rng = np.random.default_rng(0)
    dphase = np.mod(rng.normal(0, 0.3, (200, n)) + np.pi, 2 * np.pi) - np.pi
    res = evaluate(geom, dphase)
    np.testing.assert_array_equal(res.acheck, 0.0)
    residual = dphase + 2 * np.pi * res.acheck - res.param @ geom.B1.T
    assert np.abs(residual).max() < np.pi


def test_ils_is_optimal_in_the_ambiguity_metric():
    """The integer solution must beat every neighbour in {-1,0,1}^n.

    The criterion is the one the search minimises: the squared distance to the
    float solution in the metric of the ambiguity covariance ``Qahat``.  (The
    a-posteriori variance factor is a different statistic -- it fits the
    parameters without their prior -- so the two are not comparable.)
    """
    geom = realistic_geometry(seed=3)
    n = geom.Qahat.shape[0]
    rng = np.random.default_rng(7)
    dphase = np.mod(rng.normal(0, 1.0, (4, n)) + np.pi, 2 * np.pi) - np.pi
    res = evaluate(geom, dphase)

    Qi = np.linalg.inv(geom.Qahat)

    def dist(a_vec, y):
        d = y / (-2 * np.pi) - a_vec
        return d @ Qi @ d

    for i in range(4):
        best = dist(res.acheck[i], dphase[i])
        for perturb in itertools.product([-1, 0, 1], repeat=n):
            cand = res.acheck[i] + np.asarray(perturb, dtype=float)
            if np.array_equal(cand, res.acheck[i]):
                continue
            assert dist(cand, dphase[i]) >= best - 1e-6, (i, perturb)


def test_resolve_matches_the_reference_implementation():
    geom = realistic_geometry(n=6, seed=5)
    n = 6
    dphase = np.mod(np.random.default_rng(0).normal(0, 1.0, (9, n)) + np.pi, 2 * np.pi) - np.pi

    a_all = resolve_ambiguities(geom.Qahat, dphase, use_ils=True)
    model = LambdaModel(geom.Qahat)
    a_one = np.array([model.resolve_one(row)[0] for row in dphase])
    np.testing.assert_array_equal(a_all, a_one)
    assert np.all(np.isfinite(a_all))


# ---------------------------------------------------------------------------
# functional models
# ---------------------------------------------------------------------------


def test_model_definitions_match_matlab_table():
    # values transcribed from ps_model_definitions.m
    assert model_spec(1).par_index == (1, 4)
    assert model_spec(2).par_index == (1, 2, 4)
    assert model_spec(7).par_index == (1, 2, 4, 7, 8, 9, 10)
    assert npar(13) == 8
    assert model_spec(27).par_index == (1, 2, 3, 4, 5, 6)
    assert model_spec(14).defo_index == (2, 3)
    with pytest.raises(ValueError):
        model_spec(99)


def test_functional_model_columns():
    """Spot-check the design matrix against the hand-computed basis functions."""
    n = 5
    Btemp = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    h2ph = np.array([0.1, -0.1, 0.2, -0.2, 0.0])
    std = np.array([30, 0.005, 1, 0.02, 0.01, 0.01, 0.005])
    m2ph = -4 * np.pi / 0.056
    B1, Qy2 = build_functional_model(
        2, n, h2ph, Btemp, np.zeros(n), std, m2ph=m2ph, sat_vel=7550.0
    )
    np.testing.assert_allclose(B1[:, 0], m2ph * h2ph)
    np.testing.assert_allclose(B1[:, 1], m2ph)
    np.testing.assert_allclose(B1[:, 2], m2ph * Btemp)
    np.testing.assert_allclose(Qy2, [30**2, 0.005**2, 0.02**2])

    quadratic, Qy2q = build_functional_model(
        3, n, h2ph, Btemp, np.zeros(n), std, m2ph=m2ph, sat_vel=7550.0
    )
    np.testing.assert_allclose(quadratic[:, 3], m2ph * Btemp**2)
    np.testing.assert_allclose(Qy2q[3], 0.01**2)

    periodic, _ = build_functional_model(
        5, n, h2ph, Btemp, np.zeros(n), std, m2ph=m2ph, sat_vel=7550.0
    )
    np.testing.assert_allclose(periodic[:, 3], m2ph * np.sin(2 * np.pi * Btemp))
    np.testing.assert_allclose(periodic[:, 4], m2ph * (np.cos(2 * np.pi * Btemp) - 1))


def test_breakpoint_model_holds_the_rate_after_the_event():
    n = 6
    Btemp = np.array([-1.0, -0.6, -0.2, 0.2, 0.6, 1.0])
    std = np.array([30, 0.005, 1, 0.02, 0.01, 0.01, 0.005])
    m2ph = -4 * np.pi / 0.056
    B1, _ = build_functional_model(
        14, n, np.full(n, 0.1), Btemp, np.zeros(n), std,
        m2ph=m2ph, sat_vel=7550.0, breakpoint=4,
    )
    # linear1 is held constant from the breakpoint onwards
    np.testing.assert_allclose(B1[:, 2], m2ph * np.array([-1.0, -0.6, -0.2, 0.2, 0.2, 0.2]))
    # linear2 contributes only after the breakpoint
    np.testing.assert_allclose(B1[:, 3], m2ph * np.array([0, 0, 0, 0, 0.4, 0.8]))
    with pytest.raises(ValueError):
        build_functional_model(
            14, n, np.full(n, 0.1), Btemp, np.zeros(n), std,
            m2ph=m2ph, sat_vel=7550.0,
        )
