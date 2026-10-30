"""Arc network construction and variance component estimation.

Ports ``ps_form_network.m``, ``ps_construct_redundant_network.m`` and ``ps_vce.m``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay

from ..core.estimation import build_geometry, evaluate

__all__ = [
    "ArcNetwork",
    "form_network",
    "delaunay_arcs",
    "spider_arcs",
    "independent_arcs",
    "estimate_variance_components",
]


@dataclass
class ArcNetwork:
    """Arcs between PS candidates.

    ``arcs`` holds 0-based node indices, ``dphase``/``dh2ph`` are the double
    differences along each arc (``phase[j] - phase[i]``).
    """

    arcs: np.ndarray
    dphase: np.ndarray
    dh2ph: np.ndarray
    dist: np.ndarray
    Nnodes: int

    def __len__(self) -> int:
        return self.arcs.shape[0]


def delaunay_arcs(az: np.ndarray, r: np.ndarray, az_spacing: float, r_spacing: float) -> np.ndarray:
    """Delaunay triangulation -> unique (sorted) arc list, as in the MATLAB code."""
    x = np.asarray(az, dtype=float)
    y = np.asarray(r, dtype=float) * (r_spacing / az_spacing)
    if x.size < 3:
        raise ValueError("at least three PS candidates are needed for a Delaunay network")
    tri = Delaunay(np.column_stack([x, y])).simplices
    tri = np.sort(tri, axis=1)
    edges = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    return edges


def spider_arcs(
    az: np.ndarray,
    r: np.ndarray,
    az_spacing: float,
    r_spacing: float,
    Ncon: int = 8,
    max_arc_length: float = 1000.0,
) -> np.ndarray:
    """Redundant "spider" network: every node connected to its ``Ncon`` neighbours.

    Reimplementation of ``ps_construct_redundant_network.m`` using a k-d tree; the
    original divides a full circle in ``Nparts`` sectors, here the nearest
    ``Ncon`` neighbours are selected which gives a denser, more isotropic network.
    """
    from scipy.spatial import cKDTree

    x = np.asarray(az, dtype=float) * az_spacing
    y = np.asarray(r, dtype=float) * r_spacing
    tree = cKDTree(np.column_stack([x, y]))
    k = min(int(Ncon) + 1, x.size)
    dist, idx = tree.query(np.column_stack([x, y]), k=k)
    idx = np.asarray(idx).reshape(x.size, k)
    dist = np.asarray(dist).reshape(x.size, k)

    i = np.repeat(np.arange(x.size), k)
    j = idx.ravel()
    d = dist.ravel()
    mask = (j >= 0) & (j != i) & (d <= max_arc_length)
    lo = np.minimum(i[mask], j[mask]).astype(np.int64)
    hi = np.maximum(i[mask], j[mask]).astype(np.int64)
    if lo.size == 0:
        return np.array([], dtype=np.int64)
    keys = np.unique(lo * np.int64(x.size) + hi)
    return np.column_stack((keys // np.int64(x.size), keys % np.int64(x.size)))


def form_network(
    psc_az: np.ndarray,
    psc_r: np.ndarray,
    psc_phase: np.ndarray,
    psc_h2ph: np.ndarray,
    *,
    az_spacing: float,
    r_spacing: float,
    max_arc_length: float = 1000.0,
    method: str = "delaunay",
    Ncon: int = 8,
) -> ArcNetwork:
    """Build the arc network and its double differences."""
    psc_az = np.asarray(psc_az, dtype=float)
    psc_r = np.asarray(psc_r, dtype=float)
    psc_phase = np.asarray(psc_phase, dtype=float)
    psc_h2ph = np.asarray(psc_h2ph, dtype=float)

    if method == "delaunay":
        edges = delaunay_arcs(psc_az, psc_r, az_spacing, r_spacing)
    elif method == "spider":
        edges = spider_arcs(psc_az, psc_r, az_spacing, r_spacing, Ncon, max_arc_length)
    else:
        raise ValueError(f"unknown network_method {method!r}")

    if edges.size == 0:
        raise ValueError("no arcs formed")

    i, j = edges[:, 0], edges[:, 1]
    dist = np.hypot(
        az_spacing * (psc_az[i] - psc_az[j]),
        r_spacing * (psc_r[i] - psc_r[j]),
    )
    keep = dist < max_arc_length
    edges, i, j, dist = edges[keep], i[keep], j[keep], dist[keep]
    if edges.size == 0:
        raise ValueError("all arcs are longer than max_arc_length")

    dphase = psc_phase[j] - psc_phase[i]
    dh2ph = 0.5 * (psc_h2ph[j] + psc_h2ph[i])
    return ArcNetwork(arcs=edges, dphase=dphase, dh2ph=dh2ph, dist=dist, Nnodes=psc_az.size)


def independent_arcs(arcs: np.ndarray) -> np.ndarray:
    """Keep a maximally independent arc subset (no shared node).

    Equivalent to the in-place pruning loop of ``ps_vce.m``: walk the arc list and
    drop every later arc that touches an already used node.
    """
    used = np.zeros(arcs.max() + 1, dtype=bool)
    keep = []
    for idx in range(arcs.shape[0]):
        a, b = int(arcs[idx, 0]), int(arcs[idx, 1])
        if used[a] or used[b]:
            continue
        keep.append(idx)
        used[a] = used[b] = True
    return np.asarray(keep, dtype=np.int64)


def estimate_variance_components(
    arcs: np.ndarray,
    dphase: np.ndarray,
    dh2ph: np.ndarray,
    *,
    Nifgs: int,
    Btemp: np.ndarray,
    Bdop: np.ndarray | None,
    std_param: np.ndarray,
    model: int,
    m2ph: float,
    sat_vel: float,
    breakpoint: int | None = None,
    breakpoint2: int | None = None,
    method: str = "ils",
    ens_coh_threshold: float = 0.9,
    use_independent: bool = True,
) -> np.ndarray:
    """Variance-component estimation per interferogram (``ps_vce.m``).

    Returns ``sig2_est`` with ``Nifgs + 1`` entries: the master-atmosphere
    variance (zero when the master is estimated) followed by the per-interferogram
    variances.  The variance components are derived from the redundancy of the
    temporal fit on a set of mutually independent arcs.
    """
    sel = independent_arcs(arcs) if use_independent else np.arange(arcs.shape[0])
    dph = dphase[sel]
    dh = dh2ph[sel]

    geom = build_geometry(
        model, Nifgs, dh, Btemp, Bdop, std_param, None,
        m2ph=m2ph, sat_vel=sat_vel, breakpoint=breakpoint, breakpoint2=breakpoint2,
    )
    res = evaluate(geom, dph, use_ils=(method == "ils"))
    good = np.isfinite(res.varfac)
    if method == "perio":
        good = np.isfinite(res.ens_coh) & (res.ens_coh >= ens_coh_threshold)
    if not good.any():
        raise ValueError("variance component estimation failed: no valid arc")
    acheck = res.acheck[good]

    B1 = geom.B1
    Qy1 = geom.Qy1
    invQy1 = geom.invQy1
    BtQB = B1.T @ invQy1 @ B1
    Pao = np.eye(Nifgs) - B1 @ np.linalg.solve(BtQB, B1.T @ invQy1)
    QP = invQy1 @ Pao

    y = (2 * np.pi * acheck + dph[good])  # unwrapped arc phases

    # Qy1_k is rank one: 2*e_k e_k' (or 2*ones for the master term when model==1)
    if model == 1:
        Nsig = Nifgs + 1
        QPy = y @ QP.T  # (Narcs, Nifgs)
        l = np.empty((Narcs := y.shape[0], Nsig))
        l[:, 0] = 2.0 * QPy.sum(axis=1) ** 2
        l[:, 1:] = 2.0 * QPy**2
        # N[k,j] = trace(QP M_k QP M_j)
        D = QP
        N = np.empty((Nsig, Nsig))
        col = D.sum(axis=0)
        row = D.sum(axis=1)
        N[0, 0] = 4.0 * D.sum() ** 2
        N[0, 1:] = 4.0 * row * col
        N[1:, 0] = N[0, 1:]
        N[1:, 1:] = 4.0 * (D.T * D)
    else:
        Nsig = Nifgs
        QPy = y @ QP.T
        l = 2.0 * QPy**2
        D = QP
        N = 4.0 * (D.T * D)

    # N is symmetric by construction up to round-off
    N = 0.5 * (N + N.T)
    sig2 = np.linalg.solve(N, l.T)  # per arc, then averaged
    sig2_est = sig2.mean(axis=1)

    floor = (np.pi * 10 / 180.0) ** 2
    sig2_est = np.maximum(sig2_est, floor)
    if model != 1:
        sig2_est = np.concatenate(([0.0], sig2_est))
    return sig2_est
