"""Regressions for gaps that were once confirmed and are now fixed.

``GAP-S01`` (PSP payload), ``GAP-U01`` (isolated unwrapping node) and ``GAP-S02``
(side-lobe declustering) all have fixes with evidence; these tests keep them fixed.
"""
import numpy as np
import pytest
from depsipy.steps.selection import _decluster, PSCCandidate, SelectionResult, write_psc_files
from depsipy.steps.unwrapping import unwrap_network


def test_psp_export_has_observation_payload(tmp_path):
    cand = PSCCandidate(1, 1, 1, 1, np.array([0.3]), np.array([0.2]), np.array([1.]), 0.1)
    result = SelectionResult(
        np.array([[1, 2]]), np.array([[1, 2]]), np.array([1]), np.array([1]),
        psc=[[cand]], psp_az=[np.array([2])], psp_r=[np.array([2])],
        psp_phase=[np.array([[0.4]])], psp_h2ph=[np.array([[0.5]])],
        psp_amp=[np.array([[1.0]])], psp_amp_disp=[np.array([0.1])],
    )
    write_psc_files(result, tmp_path, 'tiny', 1)
    a = np.fromfile(tmp_path/'tiny_psp_2orig_sel1.raw', dtype=np.float64).reshape(1, 6)
    assert np.any(a[:, 4:] != 0)


def test_declustering_keeps_stronger_candidate():
    result = _decluster(np.array([0, 0]), np.array([0, 1]), np.zeros((2, 4)), np.array([10., 1.]), 0.95)
    np.testing.assert_array_equal(result, [0])


def test_unwrapping_allows_isolated_node():
    arcs = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    result = unwrap_network(arcs, np.zeros((6, 2)), np.zeros((6, 2)), np.ones(6), 5)
    assert np.isnan(result.phase_unw[4]).all()
    assert np.isfinite(result.phase_unw[:4]).all()
