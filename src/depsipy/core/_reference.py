"""Reference (pure numpy) implementation of the LAMBDA routines.

These functions mirror ``ldldecom.m``, ``decorrel_freek.m``, ``ps_chistart.m``
and ``bk_lsearch_count.m`` one-to-one and are the *specification* the compiled
kernels in :mod:`depsipy.core._kernels` are tested against.  They are also used
when numba is not installed.
"""

from __future__ import annotations

import numpy as np

__all__ = ["ldldecom", "decorrel", "chistart", "lsearch"]
_PY_ONLY = True


def ldldecom(Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """LtDL decomposition of a symmetric positive definite matrix.

    Returns ``L`` (unit lower triangular) and ``D`` (diagonal) such that
    ``Q = L.T @ diag(D) @ L``.

    Note the orientation: the DePSI/LAMBDA variant eliminates the *last*
    variable first, so the factors are applied as ``L.T @ diag(D) @ L`` (this
    is ``L.'*diag(D)*L`` in ``ps_chistart.m``).
    """
    Q = np.array(Q, dtype=float, copy=True)
    n = Q.shape[0]
    if Q.shape[0] != Q.shape[1]:
        raise ValueError("Q must be square")
    L = np.zeros((n, n))
    D = np.zeros(n)
    for i in range(n - 1, -1, -1):
        D[i] = Q[i, i]
        if D[i] <= 0:
            raise ValueError("matrix is not positive definite")
        L[i, : i + 1] = Q[i, : i + 1] / np.sqrt(Q[i, i])
        for j in range(i):
            Q[j, : j + 1] -= L[i, : j + 1] * L[i, j]
        L[i, : i + 1] /= L[i, i]
    if np.any(D < 1e-10):
        raise ValueError("matrix is not positive definite")
    return L, D


def decorrel(Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """LAMBDA decorrelation; returns ``(Z, L, D)``.

    ``Z`` is integer unimodular, ``L``/``D`` the LtDL factors of the decorrelated
    matrix ``Z.T @ Q @ Z``.
    """
    Q = np.asarray(Q, dtype=float)
    n = Q.shape[0]
    L, D = ldldecom(Q)
    Zti = np.eye(n)
    i1 = n - 1  # 1-based, as in decorrel_freek.m
    sw = True
    while sw:
        ii = n  # 1-based loop variable: runs for n-1, n-2, ..., 1
        sw = False
        while (not sw) and ii > 1:
            ii -= 1
            if ii <= i1:
                for jj in range(ii + 1, n + 1):
                    mu = round(L[jj - 1, ii - 1])
                    if mu != 0:
                        L[jj - 1 :, ii - 1] = L[jj - 1 :, ii - 1] - mu * L[jj - 1 :, jj - 1]
                        Zti[:, jj - 1] = Zti[:, jj - 1] + mu * Zti[:, ii - 1]
            delta = D[ii - 1] + L[ii, ii - 1] ** 2 * D[ii]
            if delta < D[ii]:
                lam = D[ii] * L[ii, ii - 1] / delta
                eta = D[ii - 1] / delta
                D[ii - 1] = eta * D[ii]
                D[ii] = delta
                help_ = L[ii, : ii - 1] - L[ii, ii - 1] * L[ii - 1, : ii - 1]
                L[ii, : ii - 1] = lam * L[ii, : ii - 1] + eta * L[ii - 1, : ii - 1]
                L[ii - 1, : ii - 1] = help_
                L[ii, ii - 1] = lam
                help_ = L[ii + 1 :, ii - 1].copy()
                L[ii + 1 :, ii - 1] = L[ii + 1 :, ii]
                L[ii + 1 :, ii] = help_
                help_ = Zti[:, ii - 1].copy()
                Zti[:, ii - 1] = Zti[:, ii]
                Zti[:, ii] = help_
                i1 = ii
                sw = True
    Z = np.round(np.linalg.inv(Zti.T)).astype(np.int64)
    # the factors of the decorrelated matrix
    Ld, Dd = ldldecom(Z.T @ Q @ Z)
    return Z, Ld, Dd


def chistart(L: np.ndarray, D: np.ndarray, a: np.ndarray, ncands: int = 1):
    """Initial size of the search ellipsoid, as ``ps_chistart.m`` with ncands=1.

    Returns ``(chi2, afixed)`` where ``afixed`` is the bootstrapped integer
    candidate (the best of the ``n+1`` partially conditioned roundings).
    """
    n = a.size
    if ncands > n + 1:  # pragma: no cover - not used by DePSI
        raise NotImplementedError("ncands > n+1 uses the volume approximation")
    invLtDL = np.linalg.inv(L.T @ np.diag(D) @ L)
    best_chi = np.inf
    best = None
    chi = np.empty(n + 1)
    cands = np.empty((n + 1, n))
    for kk in range(n + 1):  # kk = MATLAB k = n, n-1, ..., 0  -> forced level n-kk
        forced = n - kk
        afloat = a.astype(float).copy()
        afixed = a.astype(float).copy()
        for i in range(n - 1, -1, -1):
            dw = L[i:, i] @ (afloat[i:] - afixed[i:])
            afloat[i] = afloat[i] - dw
            if i != forced:
                afixed[i] = np.round(afloat[i])
            else:
                if afloat[i] == afixed[i]:
                    afixed[i] = afixed[i] + 1.0
                else:
                    afixed[i] = np.round(afloat[i] + np.sign(afloat[i] - afixed[i]))
        q = a - afixed
        chi[kk] = q @ invLtDL @ q
        cands[kk] = afixed
    idx = int(np.argmin(chi))
    best_chi = chi[idx] + 1e-6
    best = cands[idx]
    return best_chi, best


def lsearch(
    afloat: np.ndarray,
    Linv: np.ndarray,
    Dinv: np.ndarray,
    chi2: float,
    max_backtrack: int | None = None,
):
    """Branch-and-bound integer least squares (``bk_lsearch_count.m``).

    Returns ``(afixed, sqnorm, ierr, break_flag)``; ``break_flag`` is ``True``
    when the search exceeded ``max_backtrack`` backtracks, in which case the
    caller falls back to the bootstrapped candidate.
    """
    n = afloat.size
    if max_backtrack is None:
        max_backtrack = n**3

    right = np.zeros(n + 1)
    left = np.zeros(n + 1)
    right[n] = chi2
    dq = np.empty(n)
    dq[: n - 1] = Dinv[1:n] / Dinv[: n - 1]
    dq[n - 1] = 1.0 / Dinv[n - 1]

    dist = np.zeros(n)
    endd = np.zeros(n)
    lef = np.zeros(n)

    endsearch = False
    break_flag = False
    cand_n = False
    ncan = 0
    afixed = np.zeros((n, 1))
    sqnorm = np.zeros(1)
    ierr = 0

    i = n  # 0-based "one level above the last"
    iold = i
    backtrack_count = 0

    while (not endsearch) and (not break_flag):
        i -= 1
        if iold <= i:
            lef[i] = lef[i] + Linv[i + 1, i]
        else:
            if i + 1 <= n - 1:
                lef[i] = Linv[i + 1 : n, i] @ dist[i + 1 : n]
        iold = i

        right[i] = (right[i + 1] - left[i + 1]) * dq[i]
        if right[i] < 0:
            right[i] = 0.0
        reach = np.sqrt(right[i])
        delta = afloat[i] - reach - lef[i]
        dist[i] = np.ceil(delta) - afloat[i]

        if dist[i] > reach - lef[i]:
            backtrack_count += 1
            if backtrack_count == max_backtrack:
                break_flag = True
            cand_n = False
            c_stop = False
            while (not c_stop) and (i < n - 1):
                i += 1
                if dist[i] < endd[i]:
                    dist[i] += 1
                    left[i] = (dist[i] + lef[i]) ** 2
                    c_stop = True
                    if i == n - 1:
                        cand_n = True
            if (i == n - 1) and (not cand_n):
                endsearch = True
        else:
            endd[i] = reach - lef[i] - 1
            left[i] = (dist[i] + lef[i]) ** 2

        if i == 0:
            t = chi2 - (right[0] - left[0]) * Dinv[0]
            endd[0] += 1
            while dist[0] <= endd[0]:
                if ncan < 1:
                    ncan += 1
                    afixed[:, 0] = dist + afloat
                    sqnorm[0] = t
                else:
                    if t < sqnorm[0]:
                        afixed[:, 0] = dist + afloat
                        sqnorm[0] = t
                t += (2 * (dist[0] + lef[0]) + 1) * Dinv[0]
                dist[0] += 1

            cand_n = False
            c_stop = False
            backtrack_count += 1
            if backtrack_count == max_backtrack:
                break_flag = True
            while (not c_stop) and (i < n - 1):
                i += 1
                if dist[i] < endd[i]:
                    dist[i] += 1
                    left[i] = (dist[i] + lef[i]) ** 2
                    c_stop = True
                    if i == n - 1:
                        cand_n = True
            if (i == n - 1) and (not cand_n):
                endsearch = True

    if ncan < 1:
        ierr = 1
    if break_flag:
        afixed[:, 0] = np.nan
    return afixed[:, 0], float(sqnorm[0]), ierr, break_flag


