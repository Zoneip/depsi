"""Densification (step6) -- vectorised core.

This is the W2 implementation of the planned densification step.  It keeps the
depsipy data model (numpy arrays, 1-based pixel coordinates) and reuses the
batch ILS estimator.  The 2of3/3of3 consistency check is separated so it can be
tested without running the full network.

Reference: MATLAB ``ps_densification.m``; design: ``specs/SPEC-003``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..core.estimation import estimate_arcs
from ..core.models import NPAR_MAX

__all__ = [
    "DensificationResult",
    "check_consistency",
    "nearest_connections",
    "densify_points",
]


@dataclass
class DensificationResult:
    """Accepted densification points and the evidence used to accept them."""

    az: np.ndarray  # (N,) 1-based
    r: np.ndarray
    phase: np.ndarray  # (N, Nifgs)
    h2ph: np.ndarray
    amp_disp: np.ndarray
    acheck: np.ndarray  # (N, Nifgs) densification ambiguities
    param: np.ndarray  # (N, NPAR_MAX)
    varfac: np.ndarray
    ens_coh: np.ndarray
    accepted: np.ndarray  # (N,) bool
    votes: np.ndarray  # (N,) number of agreeing references
    ref_idx: np.ndarray  # (N, Nest) network indices used per candidate


def check_consistency(
    arc_acheck: np.ndarray,
    network_acheck: np.ndarray,
    ref_idx: np.ndarray,
    mode: str = "3of3",
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the 2of3 / 3of3 rule per densification point.

    Parameters
    ----------
    arc_acheck : (N, K, Nifgs)
        Integer ambiguities of each candidate-to-reference arc.
    network_acheck : (Nnet, Nifgs)
        Integer ambiguities of the network points.
    ref_idx : (N, K)
        Network indices used by each candidate.
    mode : ``'2of3'``, ``'3of3'`` or ``'nocheck'``

    Returns
    -``(accepted, votes)``.  ``votes`` counts the references (including the
    first one) whose ambiguity agrees after transferring to the first frame.
    Invalid (NaN) arcs or negative indices do not vote.
    """
    arc = np.asarray(arc_acheck, dtype=float)
    net = np.asarray(network_acheck, dtype=float)
    refs = np.asarray(ref_idx, dtype=np.int64)
    if arc.ndim != 3 or refs.shape != arc.shape[:2]:
        raise ValueError("arc_acheck and ref_idx shapes do not match")
    n_cand, k, n_ifg = arc.shape
    if k == 0:
        return np.zeros(n_cand, dtype=bool), np.zeros(n_cand, dtype=int)

    xcheck = np.empty_like(arc)
    xcheck[:, 0, :] = arc[:, 0, :]
    for w in range(1, k):
        xcheck[:, w, :] = arc[:, w, :] + net[refs[:, w], :] - net[refs[:, 0], :]

    valid = np.isfinite(xcheck).all(axis=2)
    valid &= refs >= 0
    votes = np.ones(n_cand, dtype=int)
    votes[~valid[:, 0]] = 0
    for w in range(1, k):
        diff = np.abs(xcheck[:, w, :] - xcheck[:, 0, :]).sum(axis=1)
        votes += (valid[:, w] & (diff == 0)).astype(int)

    if mode == "2of3":
        required = 2
    elif mode == "3of3":
        required = 3
    elif mode == "nocheck":
        required = 1
    else:
        raise ValueError(f"unknown dens_check {mode!r}")
    threshold = required if k >= required else k + 1
    return votes >= threshold, votes


def nearest_connections(
    cand_az: np.ndarray,
    cand_r: np.ndarray,
    net_az: np.ndarray,
    net_r: np.ndarray,
    nest: int,
    *,
    az_spacing: float = 1.0,
    r_spacing: float = 1.0,
    max_arc_length: float = 1000.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nearest ``nest`` network points per candidate, plus arc distances.

    Returns ``(ref_idx, dist, valid)``; ``ref_idx`` is ``(N, K)`` with ``-1``
    where no neighbour is available.
    """
    cand = np.column_stack([np.asarray(cand_az, float) * az_spacing, np.asarray(cand_r, float) * r_spacing])
    net = np.column_stack([np.asarray(net_az, float) * az_spacing, np.asarray(net_r, float) * r_spacing])
    n_cand, n_net = cand.shape[0], net.shape[0]
    k = min(int(nest), n_net)
    ref = np.full((n_cand, int(nest)), -1, dtype=np.int64)
    dist = np.full((n_cand, int(nest)), np.inf)
    if n_cand == 0 or n_net == 0 or k == 0:
        return ref, dist, ref >= 0
    tree = cKDTree(net)
    d, idx = tree.query(cand, k=k)
    d = np.asarray(d).reshape(n_cand, k)
    idx = np.asarray(idx).reshape(n_cand, k)
    ref[:, :k] = idx
    dist[:, :k] = d
    valid = (ref >= 0) & (dist <= max_arc_length)
    return ref, dist, valid


def densify_points(
    net_az: np.ndarray,
    net_r: np.ndarray,
    net_phase: np.ndarray,
    net_h2ph: np.ndarray,
    net_acheck: np.ndarray,
    cand_az: np.ndarray,
    cand_r: np.ndarray,
    cand_phase: np.ndarray,
    cand_h2ph: np.ndarray,
    cand_amp_disp: np.ndarray,
    *,
    Btemp: np.ndarray,
    Bdop: np.ndarray | None,
    std_param: np.ndarray,
    model: int,
    m2ph: float,
    sat_vel: float,
    nest: int = 3,
    dens_check: str = "3of3",
    max_arc_length: float = 1000.0,
    az_spacing: float = 1.0,
    r_spacing: float = 1.0,
    method: str = "ils",
    sig2_est: np.ndarray | None = None,
    breakpoint: int | None = None,
    breakpoint2: int | None = None,
) -> DensificationResult:
    """Estimate ambiguities of candidate points against the nearest network points."""
    cand_az = np.asarray(cand_az)
    cand_r = np.asarray(cand_r)
    cand_phase = np.atleast_2d(np.asarray(cand_phase, dtype=float))
    cand_h2ph = np.atleast_2d(np.asarray(cand_h2ph, dtype=float))
    net_phase = np.atleast_2d(np.asarray(net_phase, dtype=float))
    net_h2ph = np.atleast_2d(np.asarray(net_h2ph, dtype=float))
    net_acheck = np.atleast_2d(np.asarray(net_acheck, dtype=float))
    n_cand = cand_az.size
    n_ifg = cand_phase.shape[1]
    if nest < 1:
        raise ValueError("nest must be >= 1")
    if n_cand == 0:
        return DensificationResult(
            az=np.asarray(cand_az, int), r=np.asarray(cand_r, int),
            phase=cand_phase, h2ph=cand_h2ph,
            amp_disp=np.asarray(cand_amp_disp, float),
            acheck=np.empty((0, n_ifg)), param=np.empty((0, NPAR_MAX)),
            varfac=np.empty(0), ens_coh=np.empty(0),
            accepted=np.zeros(0, dtype=bool), votes=np.zeros(0, dtype=int),
            ref_idx=np.empty((0, int(nest)), dtype=np.int64),
        )
    if np.asarray(net_az).size == 0:
        return DensificationResult(
            az=np.asarray(cand_az, int), r=np.asarray(cand_r, int),
            phase=cand_phase, h2ph=cand_h2ph,
            amp_disp=np.asarray(cand_amp_disp, float),
            acheck=np.full((n_cand, n_ifg), np.nan),
            param=np.full((n_cand, NPAR_MAX), np.nan),
            varfac=np.full(n_cand, np.nan), ens_coh=np.full(n_cand, np.nan),
            accepted=np.zeros(n_cand, dtype=bool), votes=np.zeros(n_cand, dtype=int),
            ref_idx=np.full((n_cand, int(nest)), -1, dtype=np.int64),
        )

    ref, dist, valid = nearest_connections(
        cand_az, cand_r, net_az, net_r, nest,
        az_spacing=az_spacing, r_spacing=r_spacing, max_arc_length=max_arc_length,
    )
    k = ref.shape[1]
    i_cand = np.repeat(np.arange(n_cand), k)
    i_ref = ref.ravel()
    keep = (i_ref >= 0) & valid.ravel()
    i_cand, i_ref = i_cand[keep], i_ref[keep]
    if i_cand.size == 0:
        return DensificationResult(
            az=np.asarray(cand_az, int), r=np.asarray(cand_r, int),
            phase=cand_phase, h2ph=cand_h2ph,
            amp_disp=np.asarray(cand_amp_disp, float),
            acheck=np.full((n_cand, n_ifg), np.nan),
            param=np.full((n_cand, NPAR_MAX), np.nan),
            varfac=np.full(n_cand, np.nan), ens_coh=np.full(n_cand, np.nan),
            accepted=np.zeros(n_cand, dtype=bool), votes=np.zeros(n_cand, dtype=int),
            ref_idx=ref,
        )

    dphase = cand_phase[i_cand] - net_phase[i_ref]
    dh2ph = 0.5 * (cand_h2ph[i_cand] + net_h2ph[i_ref])
    est = estimate_arcs(
        dphase,
        n_ifg,
        dh2ph,
        Btemp,
        Bdop,
        std_param,
        sig2_est,
        [int(model)],
        m2ph=m2ph,
        sat_vel=sat_vel,
        method=method,
        breakpoint=breakpoint,
        breakpoint2=breakpoint2,
    )

    arc_acheck = np.full((n_cand, k, n_ifg), np.nan)
    arc_param = np.full((n_cand, k, est.param.shape[1]), np.nan)
    arc_varfac = np.full((n_cand, k), np.nan)
    arc_coh = np.full((n_cand, k), np.nan)
    # map kept arcs back to (candidate, slot); kept arcs preserve candidate order
    slots = np.zeros(n_cand, dtype=int)
    for pos, ci in enumerate(i_cand):
        slot = slots[ci]
        slots[ci] += 1
        arc_acheck[ci, slot] = est.acheck[pos]
        arc_param[ci, slot] = est.param[pos]
        arc_varfac[ci, slot] = est.varfac[pos]
        arc_coh[ci, slot] = est.ens_coh[pos]

    accepted, votes = check_consistency(arc_acheck, net_acheck, ref, dens_check)
    acheck = arc_acheck[:, 0, :] + net_acheck[ref[:, 0], :]
    acheck[ref[:, 0] < 0] = np.nan

    return DensificationResult(
        az=np.asarray(cand_az, int),
        r=np.asarray(cand_r, int),
        phase=cand_phase,
        h2ph=cand_h2ph,
        amp_disp=np.asarray(cand_amp_disp, float),
        acheck=acheck,
        param=arc_param[:, 0, :],
        varfac=arc_varfac[:, 0],
        ens_coh=arc_coh[:, 0],
        accepted=accepted,
        votes=votes,
        ref_idx=ref,
    )
