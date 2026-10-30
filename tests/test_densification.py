"""Densification core tests (F1/W2)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from depsipy.core.models import NPAR_MAX
from depsipy.steps import densification as dens


def test_nearest_connections_match_bruteforce():
    rng = np.random.default_rng(4)
    net = rng.uniform(0, 100, (20, 2))
    cand = rng.uniform(0, 100, (7, 2))
    ref, dist, valid = dens.nearest_connections(cand[:, 0], cand[:, 1], net[:, 0], net[:, 1], 3)
    assert ref.shape == (7, 3)
    for i in range(cand.shape[0]):
        brute = np.sqrt(((net - cand[i]) ** 2).sum(axis=1))
        want = np.sort(brute)[:3]
        np.testing.assert_allclose(np.sort(dist[i]), want)
        assert ref[i, 0] == int(np.argmin(brute))
        assert valid[i].all()


def test_check_consistency_2of3_and_3of3():
    net_acheck = np.zeros((4, 2))
    refs = np.array([[0, 1, 2], [1, 2, 3], [0, 2, 3]])
    arcs = np.array([
        [[0, 0], [0, 0], [0, 0]],
        [[0, 0], [0, 0], [1, 1]],
        [[0, 0], [1, 1], [2, 2]],
    ], dtype=float)
    accepted2, votes2 = dens.check_consistency(arcs, net_acheck, refs, "2of3")
    accepted3, votes3 = dens.check_consistency(arcs, net_acheck, refs, "3of3")
    np.testing.assert_array_equal(votes2, [3, 2, 1])
    np.testing.assert_array_equal(votes3, [3, 2, 1])
    np.testing.assert_array_equal(accepted2, [True, True, False])
    np.testing.assert_array_equal(accepted3, [True, False, False])
    accepted_none, _ = dens.check_consistency(arcs[:, :1, :], net_acheck, refs[:, :1], "3of3")
    assert not accepted_none.any()


def test_check_consistency_rejects_invalid_arcs():
    net_acheck = np.zeros((3, 1))
    refs = np.array([[0, 1, 2], [-1, -1, -1]])
    arcs = np.array([[[0.], [0.], [0.]], [[np.nan], [np.nan], [np.nan]]])
    accepted, votes = dens.check_consistency(arcs, net_acheck, refs, "2of3")
    np.testing.assert_array_equal(accepted, [True, False])
    np.testing.assert_array_equal(votes, [3, 0])


def test_densify_points_uses_batch_estimator(monkeypatch):
    n_ifg = 2
    k = 3
    net_az = np.array([1, 1, 1, 1])
    net_r = np.array([1, 10, 20, 30])
    net_phase = np.zeros((4, n_ifg))
    net_h2ph = np.zeros((4, n_ifg))
    net_acheck = np.zeros((4, n_ifg))
    cand_az = np.array([1, 1, 1])
    cand_r = np.array([1, 2, 3])
    cand_phase = np.zeros((3, n_ifg))
    cand_h2ph = np.zeros((3, n_ifg))
    amp = np.array([0.1, 0.2, 0.3])

    n_arcs = 3 * k
    fake_acheck = np.zeros((n_arcs, n_ifg))
    fake_acheck[3:6, 0] = [0, 0, 1]  # candidate 1: one disagreeing reference
    fake_acheck[6:9, 0] = [0, 1, 2]  # candidate 2: only ref 0 agrees
    fake = SimpleNamespace(
        acheck=fake_acheck,
        param=np.zeros((n_arcs, NPAR_MAX)),
        varfac=np.ones(n_arcs),
        ens_coh=np.full(n_arcs, 0.9),
    )
    monkeypatch.setattr(dens, "estimate_arcs", lambda *a, **kw: fake)

    out = dens.densify_points(
        net_az, net_r, net_phase, net_h2ph, net_acheck,
        cand_az, cand_r, cand_phase, cand_h2ph, amp,
        Btemp=np.array([0.0, 1.0]), Bdop=None, std_param=np.ones(7),
        model=1, m2ph=1.0, sat_vel=1.0, nest=k, dens_check="2of3",
        az_spacing=1.0, r_spacing=1.0, max_arc_length=1000.0,
    )
    np.testing.assert_array_equal(out.accepted, [True, True, False])
    np.testing.assert_array_equal(out.votes, [3, 2, 1])
    assert out.acheck.shape == (3, n_ifg)
    assert out.param.shape == (3, NPAR_MAX)

    out3 = dens.densify_points(
        net_az, net_r, net_phase, net_h2ph, net_acheck,
        cand_az, cand_r, cand_phase, cand_h2ph, amp,
        Btemp=np.array([0.0, 1.0]), Bdop=None, std_param=np.ones(7),
        model=1, m2ph=1.0, sat_vel=1.0, nest=k, dens_check="3of3",
        az_spacing=1.0, r_spacing=1.0, max_arc_length=1000.0,
    )
    np.testing.assert_array_equal(out3.accepted, [True, False, False])


def test_densify_empty_inputs():
    zeros = np.zeros((0, 2))
    out = dens.densify_points(
        np.zeros(0), np.zeros(0), zeros, zeros, zeros,
        np.zeros(0), np.zeros(0), zeros, zeros, np.zeros(0),
        Btemp=np.array([0.0, 1.0]), Bdop=None, std_param=np.ones(7),
        model=1, m2ph=1.0, sat_vel=1.0, nest=3,
    )
    assert out.accepted.size == 0
    assert out.acheck.shape == (0, 2)


def test_densify_without_network_rejects_everything():
    n_ifg = 2
    out = dens.densify_points(
        np.zeros(0), np.zeros(0), np.zeros((0, n_ifg)), np.zeros((0, n_ifg)),
        np.zeros((0, n_ifg)),
        np.array([1, 2]), np.array([1, 2]), np.zeros((2, n_ifg)),
        np.zeros((2, n_ifg)), np.array([0.1, 0.2]),
        Btemp=np.array([0.0, 1.0]), Bdop=None, std_param=np.ones(7),
        model=1, m2ph=1.0, sat_vel=1.0, nest=3,
    )
    assert not out.accepted.any()
