"""Parity tests for the vectorised loops (W1/W2/W3).

Each optimised routine is compared against the explicit loop implementation it
replaced.  The reference code is kept here, not in production, so the parity
contract is visible and cannot silently drift.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree

from depsipy.steps.network import spider_arcs
from depsipy.steps.selection import _decluster
from depsipy.steps.unwrapping import (
    _normal_matrix,
    _resolve_integer_field_batched,
    _resolve_integer_field_loop,
    _solve_multi,
)


def _ref_spider_arcs(az, r, az_spacing, r_spacing, Ncon, max_arc_length):
    x = np.asarray(az, dtype=float) * az_spacing
    y = np.asarray(r, dtype=float) * r_spacing
    tree = cKDTree(np.column_stack([x, y]))
    k = min(int(Ncon) + 1, x.size)
    dist, idx = tree.query(np.column_stack([x, y]), k=k)
    idx = np.asarray(idx).reshape(x.size, k)
    dist = np.asarray(dist).reshape(x.size, k)
    edges = set()
    for i in range(x.size):
        for j, d in zip(idx[i], dist[i]):
            if j == i or d > max_arc_length:
                continue
            edges.add((min(i, int(j)), max(i, int(j))))
    if not edges:
        return np.array([], dtype=np.int64)
    return np.array(sorted(edges), dtype=np.int64)


@pytest.mark.parametrize("n, ncon", [(1, 8), (20, 0), (20, 8), (50, 4)])
def test_spider_arcs_matches_loop(n, ncon):
    rng = np.random.default_rng(10 + n + ncon)
    az = rng.uniform(0, 100, n)
    r = rng.uniform(0, 100, n)
    got = spider_arcs(az, r, 1.0, 1.0, Ncon=ncon, max_arc_length=80.0)
    ref = _ref_spider_arcs(az, r, 1.0, 1.0, ncon, 80.0)
    np.testing.assert_array_equal(got, ref)


def test_spider_arcs_duplicate_coordinates():
    az = np.array([0.0, 0.0, 1.0, 1.0])
    r = np.array([0.0, 0.0, 1.0, 1.0])
    got = spider_arcs(az, r, 1.0, 1.0, Ncon=3, max_arc_length=10.0)
    ref = _ref_spider_arcs(az, r, 1.0, 1.0, 3, 10.0)
    np.testing.assert_array_equal(got, ref)


def _small_network(n_nodes=16, n_ifgs=5, seed=3):
    """A connected grid-like network with a reference node removed."""
    rng = np.random.default_rng(seed)
    side = 4
    pairs = []
    for a in range(side):
        for b in range(side):
            node = a * side + b
            if b + 1 < side:
                pairs.append((node, node + 1))
            if a + 1 < side:
                pairs.append((node, node + side))
    arcs = np.array(pairs, dtype=np.int64)
    weights = rng.uniform(0.5, 2.0, arcs.shape[0])
    free = np.ones(n_nodes, dtype=bool)
    free[0] = False
    Af, N, AtW = _normal_matrix(arcs, n_nodes, weights, free)
    solve = spla.factorized(N.tocsc())
    Y = rng.integers(-2, 3, (arcs.shape[0], n_ifgs)).astype(float)
    # Perturb a few integers so the correction loop has work to do.
    Y[rng.choice(arcs.shape[0], size=max(1, arcs.shape[0] // 8), replace=False), 0] += 1
    return Y, Af, AtW, solve, weights


@pytest.mark.parametrize("max_fixups", [3, 10, 200])
def test_batched_integer_correction_matches_loop(max_fixups):
    Y, Af, AtW, solve, weights = _small_network()
    Y_ref, n_ref = _resolve_integer_field_loop(
        Y, Af, AtW, solve, weights, max_fixups, 1e-10
    )
    Y_vec, n_vec = _resolve_integer_field_batched(
        Y, Af, AtW, solve, weights, max_fixups, 1e-10
    )
    np.testing.assert_allclose(Y_vec, Y_ref, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(n_vec, n_ref)


def test_solve_multi_matches_single_column_and_chunks():
    rng = np.random.default_rng(11)
    AtW = sp.random(30, 12, density=0.3, format="csr", random_state=rng)
    M = rng.normal(size=(12, 5))
    solve = lambda b: np.asarray(b)
    direct = _solve_multi(solve, AtW, M)
    chunked = _solve_multi(solve, AtW, M, max_bytes=1)
    np.testing.assert_allclose(chunked, direct)
    single = _solve_multi(solve, AtW, M[:, 0])
    np.testing.assert_allclose(single, direct[:, 0])
    empty = _solve_multi(solve, AtW, np.empty((12, 0)))
    assert empty.shape == (30, 0)


def test_decluster_keeps_strongest_drops_coherent_weaker():
    # same row: the strongest is the reference, its coherent neighbour is dropped
    rows = np.array([0, 0])
    cols = np.array([0, 1])
    theta = np.zeros((2, 4))
    got = _decluster(rows, cols, theta, np.array([10.0, 1.0]), 0.95)
    np.testing.assert_array_equal(got, [0])


def test_decluster_keeps_incoherent_candidates():
    # middle candidate is coherent with the strongest, the third is not
    rows = np.array([0, 0, 0])
    cols = np.array([0, 1, 2])
    theta = np.array([[0.0] * 4, [0.0] * 4, [0.0, np.pi, 0.0, np.pi]])
    got = _decluster(rows, cols, theta, np.array([10.0, 2.0, 1.0]), 0.95)
    np.testing.assert_array_equal(got, [0, 2])


def test_decluster_azimuth_sweep_and_duplicates():
    # same column, coherent: only the strongest survives
    rows = np.array([0, 1, 2])
    cols = np.array([0, 0, 0])
    theta = np.zeros((3, 3))
    got = _decluster(rows, cols, theta, np.array([1.0, 5.0, 2.0]), 0.95)
    np.testing.assert_array_equal(got, [1])

    # duplicated coordinates must not be "nearest neighbours" of each other
    rows = np.array([0, 0, 0, 1])
    cols = np.array([0, 0, 1, 0])
    theta = np.zeros((4, 3))
    amp = np.array([1.0, 2.0, 3.0, 4.0])
    out = _decluster(rows, cols, theta, amp, 0.95)
    assert out.size >= 1 and np.issubdtype(out.dtype, np.integer)
