"""Spatial unwrapping of the PS network (port of ``ps_spatial_unwrapping.m``).

The DePSI scheme is a two-step network adjustment:

1. **Ambiguity adjustment.**  The integer ambiguities estimated per arc are
   inconsistent across the network because of estimation errors.  A weighted
   least-squares adjustment of the *integer* field is solved (``acheck`` per PS)
   and individual integers are corrected one at a time (largest residual first)
   until the overall model test passes.
2. **Phase adjustment.**  With consistent integers, the same normal equations are
   solved for the *wrapped* double differences, which yields the unwrapped phase
   at every PS relative to the reference point of its connected group.

Both steps use one sparse normal matrix, so this implementation factors it once
and solves all interferograms as multiple right-hand sides -- the MATLAB code
re-factors inside a per-arc loop.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse import csgraph

__all__ = [
    "UnwrapResult",
    "connected_groups",
    "remove_untestable",
    "unwrap_network",
]


@dataclass
class UnwrapResult:
    acheck: np.ndarray  # (Npsc, Nifgs) integer ambiguities, NaN where excluded
    phase_unw: np.ndarray  # (Npsc, Nifgs) unwrapped phase [rad], NaN where excluded
    ref_nodes: np.ndarray  # (Nref,) 0-based node indices used as references
    ref_coords: np.ndarray  # (Nref, 2) (az, r) of the references
    groups: np.ndarray  # (Npsc,) group id, 0 = not connected
    arcs: np.ndarray  # (Narcs, 2) kept arcs
    removed: np.ndarray  # (Nremoved, 2) removed arcs
    adapted: int  # number of corrected integers
    n_fixups: np.ndarray  # (Nifgs,) corrections applied per interferogram

    @property
    def n_ref(self) -> int:
        return int(self.ref_nodes.size)


def connected_groups(arcs: np.ndarray, Nnodes: int) -> np.ndarray:
    """Connected component id per node (0 = isolated, 1.. = groups)."""
    if arcs.size == 0:
        return np.zeros(Nnodes, dtype=np.int64)
    A = sp.coo_matrix(
        (np.ones(arcs.shape[0]), (arcs[:, 0], arcs[:, 1])), shape=(Nnodes, Nnodes)
    )
    ncomp, labels = csgraph.connected_components(A, directed=False)
    # nodes without any arc keep label 0 -> shift so that groups start at 1
    has_arc = np.zeros(Nnodes, dtype=bool)
    has_arc[arcs.ravel()] = True
    groups = np.where(has_arc, labels + 1, 0)
    return groups


def remove_untestable(arcs: np.ndarray, min_connections: int = 3):
    """Iteratively drop arcs of nodes with fewer than ``min_connections`` arcs.

    Returns ``(kept_arcs, removed_arcs)``; a node needs at least three arcs to be
    testable in a network adjustment (two would make it a dangling chain).
    """
    arcs = np.asarray(arcs, dtype=np.int64)
    removed = []
    while arcs.shape[0]:
        counts = np.bincount(arcs.ravel())
        bad = np.nonzero(counts < min_connections)[0]
        if bad.size == 0:
            break
        mask = np.isin(arcs[:, 0], bad) | np.isin(arcs[:, 1], bad)
        if not mask.any():
            break
        removed.append(arcs[mask])
        arcs = arcs[~mask]
    removed_all = np.vstack(removed) if removed else np.zeros((0, 2), dtype=np.int64)
    return arcs, removed_all


def _normal_matrix(arcs, Nnodes, weights, free_cols):
    """Build ``A' W A`` and ``A' W`` restricted to the free columns.

    ``A`` has -1 in the first node and +1 in the second node of each arc, so the
    adjustment estimates the PS phase differences.  ``free_cols`` is a boolean
    mask over the nodes (references and nodes without arcs are excluded); their
    columns are dropped from the unknowns, which is how the references are moved
    to the right-hand side.  A node without arcs must be dropped as well, or its
    zero column makes the normal matrix singular.
    """
    narcs = arcs.shape[0]
    rows = np.repeat(np.arange(narcs), 2)
    cols = arcs.ravel()
    vals = np.tile(np.array([-1.0, 1.0]), narcs)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(narcs, Nnodes))

    Af = A[:, np.asarray(free_cols, dtype=bool)]
    W = sp.diags(weights) if np.ndim(weights) else sp.identity(narcs) * weights
    N = (Af.T @ W @ Af).tocsc()
    AtW = (Af.T @ W).tocsr()
    return Af, N, AtW


def _resolve_integer_field_loop(Y, Af, AtW, solve, weights, max_fixups, tol):
    """Reference per-IFG correction loop, retained for parity tests."""
    Y = np.asarray(Y, dtype=float).copy()
    n_fixups = np.zeros(Y.shape[1], dtype=int)
    for v in range(Y.shape[1]):
        y = Y[:, v].copy()
        skip = 0
        prev = -1
        for _ in range(max_fixups):
            rhs = AtW @ y
            x = solve(rhs)
            e = y - Af @ x
            omt = float(e @ (weights * e))
            if omt <= tol:
                break
            if skip == 0:
                idx = int(np.argmax(np.abs(e)))
            else:
                order = np.argsort(np.abs(e))
                idx = int(order[-min(skip + 1, order.size)])
            if idx == prev:
                skip += 1
                continue
            skip = 0
            step = -np.round(e[idx])
            if step == 0:
                step = -np.sign(e[idx]) if e[idx] != 0 else 1.0
            y[idx] += step
            n_fixups[v] += 1
            prev = idx
        Y[:, v] = y
    return Y, n_fixups


def _resolve_integer_field_batched(Y, Af, AtW, solve, weights, max_fixups, tol):
    """Correction loop with all interferograms advanced in one right-hand side.

    Each IFG keeps its own sequential skip/prev state; only the outer IFG loop
    is batched.  This preserves the per-IFG decision rules while turning
    ``Nifgs`` sparse solves per iteration into one.
    """
    Y = np.asarray(Y, dtype=float).copy()
    n_arcs, n_ifgs = Y.shape
    n_fixups = np.zeros(n_ifgs, dtype=int)
    if n_ifgs == 0 or max_fixups <= 0:
        return Y, n_fixups

    active = np.arange(n_ifgs)
    skip = np.zeros(n_ifgs, dtype=int)
    prev = np.full(n_ifgs, -1, dtype=int)
    iters = np.zeros(n_ifgs, dtype=int)

    while active.size:
        cols = active
        rhs = AtW @ Y[:, cols]
        x = solve(rhs)
        e = Y[:, cols] - Af @ x
        omt = np.sum(e * (weights[:, None] * e), axis=0)
        converged = omt <= tol
        iters += 1

        idx = np.argmax(np.abs(e), axis=0)
        if np.any(skip > 0):
            sub = np.nonzero(skip > 0)[0]
            order = np.argsort(np.abs(e[:, sub]), axis=0)
            kth = np.minimum(skip[sub] + 1, n_arcs)
            idx = np.asarray(idx)
            idx[sub] = order[n_arcs - kth, np.arange(sub.size)]

        same = idx == prev
        update = ~same & ~converged
        if update.any():
            pos = np.nonzero(update)[0]
            e_idx = e[idx[update], pos]
            step = -np.round(e_idx)
            step = np.where(step == 0, np.where(e_idx != 0, -np.sign(e_idx), 1.0), step)
            Y[idx[update], cols[update]] += step
            n_fixups[cols[update]] += 1
            prev[update] = idx[update]
            skip[update] = 0
        skip[same & ~converged] += 1

        done = converged | (iters >= max_fixups)
        active = cols[~done]
        skip = skip[~done]
        prev = prev[~done]
        iters = iters[~done]

    return Y, n_fixups


def _solve_multi(solve, AtW, M, max_bytes: int = 256 * 1024 * 1024):
    """Solve ``AtW @ M`` for one or more right-hand sides, chunking if needed."""
    M = np.asarray(M, dtype=float)
    if M.ndim == 1:
        return solve(AtW @ M)
    n_free = AtW.shape[0]
    n_cols = M.shape[1]
    if n_cols == 0:
        return np.zeros((n_free, 0), dtype=float)
    if n_free * n_cols * 8 <= max_bytes:
        return solve(AtW @ M)
    out = np.empty((n_free, n_cols), dtype=float)
    chunk = max(1, max_bytes // max(1, n_free * 8))
    for start in range(0, n_cols, chunk):
        stop = min(start + chunk, n_cols)
        out[:, start:stop] = solve(AtW @ M[:, start:stop])
    return out


def unwrap_network(
    arcs: np.ndarray,
    dphase: np.ndarray,
    dacheck: np.ndarray,
    varfac: np.ndarray,
    Nnodes: int,
    *,
    ref_cn=None,
    weighted: bool = True,
    varfac_threshold: float = 3.0,
    ens_coh: np.ndarray | None = None,
    ens_coh_threshold: float = 0.9,
    method: str = "ils",
    Nlines: int | None = None,
    Npixels: int | None = None,
    az_spacing: float = 1.0,
    r_spacing: float = 1.0,
    ps_az: np.ndarray | None = None,
    ps_r: np.ndarray | None = None,
    max_fixups: int | None = None,
    omt_tolerance: float = 1e-10,
    min_connections: int = 3,
) -> UnwrapResult:
    """Adjust the network and unwrap.

    Parameters
    ----------
    arcs : (Narcs, 2) int
        0-based node indices, ``dphase`` the wrapped double difference
        ``phase[j] - phase[i]``.
    dacheck : (Narcs, Nifgs) int
        Estimated integer ambiguities per arc.
    varfac : (Narcs,) float
        A-posteriori variance factors, used both as the arc weight
        (``w = 1/varfac**2``) and as the a-priori quality filter.
    """
    Nifgs = dphase.shape[1]
    arcs = np.asarray(arcs, dtype=np.int64)
    dphase = np.asarray(dphase, dtype=float)
    dacheck = np.asarray(dacheck, dtype=float)
    varfac = np.asarray(varfac, dtype=float)

    if method == "perio":
        if ens_coh is None:
            raise ValueError("method='perio' needs the ensemble coherence")
        keep = np.isfinite(ens_coh) & (ens_coh >= ens_coh_threshold)
    else:
        keep = np.isfinite(varfac) & (varfac <= varfac_threshold)
    removed = [arcs[~keep]]
    arcs = arcs[keep]
    dphase = dphase[keep]
    dacheck = dacheck[keep]
    varfac = varfac[keep]
    if arcs.shape[0] == 0:
        raise ValueError("no arc survived the quality filter")

    # ``remove_untestable`` returns the surviving arcs only; the parallel arrays
    # still hold every arc of this step, so filter all of them with one mask that
    # is built on the *pre-removal* arc list.
    arcs_in = arcs
    arcs, dropped = remove_untestable(arcs, min_connections)
    removed.append(dropped)
    if arcs.shape[0] == 0:
        raise ValueError("no testable arc left after removing dangling nodes")

    if dropped.size:
        if Nnodes >= 3_000_000_000:
            raise ValueError("Nnodes too large for int64 arc keys")
        dropped_key = dropped[:, 0] * np.int64(Nnodes) + dropped[:, 1]
        arc_key = arcs_in[:, 0] * np.int64(Nnodes) + arcs_in[:, 1]
        mask = ~np.isin(arc_key, dropped_key)
        arcs = arcs_in[mask]
        dphase, dacheck, varfac = dphase[mask], dacheck[mask], varfac[mask]

    groups = connected_groups(arcs, Nnodes)

    # ---- reference point per group ---------------------------------------
    is_ref = np.zeros(Nnodes, dtype=bool)
    ref_nodes = []
    if ref_cn is not None and ps_az is not None and np.all(np.isfinite(ref_cn)):
        match = np.nonzero((ps_az == ref_cn[0]) & (ps_r == ref_cn[1]))[0]
        if match.size == 0:
            raise ValueError(f"reference point {tuple(ref_cn)} is not a PS")
        ref_nodes.append(int(match[0]))
    else:
        # the node of the arc with the smallest variance factor
        best = int(np.argmin(varfac))
        ref_nodes.append(int(arcs[best, 0]))

    present = np.unique(arcs.ravel())
    for g in np.unique(groups[present]):
        if groups[ref_nodes[0]] == g:
            continue
        nodes = present[groups[present] == g]
        if nodes.size == 0:
            continue
        if ps_az is not None and Nlines is not None:
            ca = Nlines * az_spacing / 2.0
            cr = Npixels * r_spacing / 2.0
            d = (ps_az[nodes] * az_spacing - ca) ** 2 + (ps_r[nodes] * r_spacing - cr) ** 2
            ref_nodes.append(int(nodes[int(np.argmin(d))]))
        else:
            ref_nodes.append(int(nodes[0]))
    is_ref[np.asarray(ref_nodes, dtype=int)] = True

    ref_coords = (
        np.column_stack([ps_az[ref_nodes], ps_r[ref_nodes]])
        if ps_az is not None
        else np.zeros((len(ref_nodes), 2))
    )

    # ---- normal equations -------------------------------------------------
    weights = 1.0 / np.maximum(varfac, 1e-6) ** 2 if weighted else np.ones_like(varfac)
    free = ~is_ref & (groups != 0)
    Af, N, AtW = _normal_matrix(arcs, Nnodes, weights, free)
    solve = spla.factorized(N.tocsc())
    Nfree = int(free.sum())

    # ---- step 1: inconsistent integers -----------------------------------
    if max_fixups is None:
        max_fixups = arcs.shape[0]
    Y, n_fixups = _resolve_integer_field_batched(
        dacheck, Af, AtW, solve, weights, max_fixups, omt_tolerance
    )

    # ---- final integer + phase adjustment --------------------------------
    Xa = _solve_multi(solve, AtW, Y) if Nifgs else np.zeros((Nfree, 0))
    Xp = _solve_multi(solve, AtW, dphase)

    acheck = np.full((Nnodes, Nifgs), np.nan)
    phase_unw = np.full((Nnodes, Nifgs), np.nan)
    acheck[free] = np.round(Xa)
    phase_unw[free] = np.round(Xa) * 2 * np.pi + Xp
    acheck[is_ref] = 0.0
    phase_unw[is_ref] = 0.0

    return UnwrapResult(
        acheck=acheck,
        phase_unw=phase_unw,
        ref_nodes=np.asarray(ref_nodes, dtype=int),
        ref_coords=np.asarray(ref_coords),
        groups=groups,
        arcs=arcs,
        removed=np.vstack([r for r in removed if r.size]) if any(r.size for r in removed) else np.zeros((0, 2), np.int64),
        adapted=int(n_fixups.sum()),
        n_fixups=n_fixups,
    )
