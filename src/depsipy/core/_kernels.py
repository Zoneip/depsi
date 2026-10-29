"""Compiled kernels for the LAMBDA ambiguity resolution loop.

The public function :func:`resolve_all` fills an ``(Narcs, Nifgs)`` array with
the integer ambiguities of every arc.  With numba installed the branch-and-bound
search runs as machine code, parallelised over arcs; without numba it falls back
to the pure-Python reference implementation in :mod:`._reference`.

The arithmetic is deliberately kept in the same operation order as MATLAB so
that results are bit-comparable: the only differences are that the search
ellipsoid bound is clipped at zero (MATLAB would produce NaN in the ``sqrt``)
and that the FvL "max_backtrack" guard is applied per arc.
"""

from __future__ import annotations

import numpy as np

from ._reference import chistart as _py_chistart
from ._reference import lsearch as _py_lsearch

__all__ = ["HAVE_NUMBA", "resolve_all", "nb_chistart", "nb_lsearch"]

try:  # pragma: no cover - depends on the environment
    import numba as _nb
    from numba import njit, prange

    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    _nb = None
    HAVE_NUMBA = False


# ---------------------------------------------------------------------------
# numba kernels: 1:1 translations of the MATLAB loops
# ---------------------------------------------------------------------------
if HAVE_NUMBA:

    @njit(cache=True, inline="always")
    def _mround(x):
        """MATLAB round(): halves away from zero."""
        if x >= 0.0:
            return np.floor(x + 0.5)
        return np.ceil(x - 0.5)

    @njit(cache=True, inline="always")
    def _chistart_nb(L, D, a, invLtDL, out):
        """Best of the n+1 partially conditioned roundings (ncands = 1)."""
        n = a.shape[0]
        afloat = np.empty(n)
        afixed = np.empty(n)
        q = np.empty(n)
        best_chi = 0.0
        found = False
        for kk in range(n + 1):
            forced = n - kk
            for i in range(n):
                afloat[i] = a[i]
                afixed[i] = a[i]
            for i in range(n - 1, -1, -1):
                dw = 0.0
                for j in range(i, n):
                    dw += L[j, i] * (afloat[j] - afixed[j])
                afloat[i] = afloat[i] - dw
                if i != forced:
                    afixed[i] = _mround(afloat[i])
                else:
                    if afloat[i] == afixed[i]:
                        afixed[i] = afixed[i] + 1.0
                    else:
                        s = 1.0 if (afloat[i] - afixed[i]) >= 0 else -1.0
                        afixed[i] = _mround(afloat[i] + s)
            for i in range(n):
                q[i] = a[i] - afixed[i]
            chi = 0.0
            for i in range(n):
                acc = 0.0
                for j in range(n):
                    acc += invLtDL[i, j] * q[j]
                chi += q[i] * acc
            if (not found) or (chi < best_chi):
                found = True
                best_chi = chi
                for i in range(n):
                    out[i] = afixed[i]
        return best_chi + 1e-6

    @njit(cache=True)
    def _lsearch_nb(afloat, Linv, Dinv, chi2, max_backtrack, afixed):
        """Branch-and-bound search; returns (sqnorm, ierr, break_flag)."""
        n = afloat.shape[0]
        right = np.zeros(n + 1)
        left = np.zeros(n + 1)
        right[n] = chi2
        dq = np.zeros(n)
        for i in range(n - 1):
            dq[i] = Dinv[i + 1] / Dinv[i]
        dq[n - 1] = 1.0 / Dinv[n - 1]

        dist = np.zeros(n)
        endd = np.zeros(n)
        lef = np.zeros(n)

        endsearch = False
        break_flag = False
        cand_n = False
        ncan = 0
        sqnorm = 0.0
        ierr = 0
        backtrack_count = 0

        i = n
        iold = i
        while (not endsearch) and (not break_flag):
            i -= 1
            if iold <= i:
                lef[i] = lef[i] + Linv[i + 1, i]
            else:
                if i + 1 <= n - 1:
                    acc = 0.0
                    for j in range(i + 1, n):
                        acc += Linv[j, i] * dist[j]
                    lef[i] = acc
            iold = i

            right[i] = (right[i + 1] - left[i + 1]) * dq[i]
            if right[i] < 0.0:
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
                        dist[i] = dist[i] + 1
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
                endd[0] = endd[0] + 1
                while dist[0] <= endd[0]:
                    if ncan < 1:
                        ncan = 1
                        for j in range(n):
                            afixed[j] = dist[j] + afloat[j]
                        sqnorm = t
                    else:
                        if t < sqnorm:
                            for j in range(n):
                                afixed[j] = dist[j] + afloat[j]
                            sqnorm = t
                    t = t + (2.0 * (dist[0] + lef[0]) + 1.0) * Dinv[0]
                    dist[0] = dist[0] + 1

                cand_n = False
                c_stop = False
                backtrack_count += 1
                if backtrack_count == max_backtrack:
                    break_flag = True
                while (not c_stop) and (i < n - 1):
                    i += 1
                    if dist[i] < endd[i]:
                        dist[i] = dist[i] + 1
                        left[i] = (dist[i] + lef[i]) ** 2
                        c_stop = True
                        if i == n - 1:
                            cand_n = True
                if (i == n - 1) and (not cand_n):
                    endsearch = True

        if ncan < 1:
            ierr = 1
        if break_flag:
            for j in range(n):
                afixed[j] = np.nan
        return sqnorm, ierr, break_flag

    @njit(cache=True, parallel=True, nogil=True)
    def _resolve_all_nb(dphase, Zt, invLtDL, L, D, invLin, invD, invZT, use_ils, max_backtrack, out):
        npnts, n = dphase.shape
        zfloat = np.empty(n)
        zboot = np.empty(n)
        zcheck = np.empty(n)
        afixed = np.empty(n)
        for p in prange(npnts):
            for i in range(n):
                zfloat[i] = dphase[p, i] / (-2.0 * np.pi)
            # zfloat = Zt @ afloat
            for i in range(n):
                acc = 0.0
                for j in range(n):
                    acc += Zt[i, j] * zfloat[j]
                afixed[i] = acc
            for i in range(n):
                zfloat[i] = afixed[i]

            chi2 = _chistart_nb(L, D, zfloat, invLtDL, zboot)

            if use_ils:
                _lsearch_nb(zfloat, invLin, invD, chi2, max_backtrack, zcheck)
                # fall back to the bootstrapped candidate when the search bailed out
                brk = False
                for i in range(n):
                    if np.isnan(zcheck[i]):
                        brk = True
                if brk:
                    for i in range(n):
                        zcheck[i] = zboot[i]
            else:
                for i in range(n):
                    zcheck[i] = zboot[i]

            # acheck = round(invZT @ zcheck)
            for i in range(n):
                acc = 0.0
                for j in range(n):
                    acc += invZT[i, j] * zcheck[j]
                out[p, i] = _mround(acc)

    def nb_chistart(L, D, invLtDL, a):
        out = np.zeros(a.size)
        chi2 = _chistart_nb(
            np.ascontiguousarray(L), np.ascontiguousarray(D),
            np.ascontiguousarray(a.astype(float)),
            np.ascontiguousarray(np.asarray(invLtDL, dtype=float)), out,
        )
        return chi2, out

    def nb_lsearch(afloat, Linv, Dinv, chi2, max_backtrack=None):
        n = afloat.size
        if max_backtrack is None:
            max_backtrack = n**3
        afixed = np.zeros(n)
        sq, ierr, brk = _lsearch_nb(
            np.ascontiguousarray(afloat.astype(float)),
            np.ascontiguousarray(Linv), np.ascontiguousarray(Dinv),
            float(chi2), int(max_backtrack), afixed,
        )
        return afixed, sq, ierr, brk

else:  # pragma: no cover
    nb_chistart = None
    nb_lsearch = None


# ---------------------------------------------------------------------------
# dispatcher
# ---------------------------------------------------------------------------


def resolve_all(dphase, Zt, invLtDL, L, D, invL, invD, invZT, use_ils, max_backtrack, out):
    """Resolve every arc in ``dphase`` into ``out`` (integer ambiguities)."""
    dphase = np.ascontiguousarray(dphase, dtype=float)
    out_arr = np.ascontiguousarray(out, dtype=float)
    n = dphase.shape[1]

    if HAVE_NUMBA and dphase.shape[0] >= 8:
        _resolve_all_nb(
            dphase,
            np.ascontiguousarray(Zt, dtype=float),
            np.ascontiguousarray(invLtDL, dtype=float),
            np.ascontiguousarray(L, dtype=float),
            np.ascontiguousarray(D, dtype=float),
            np.ascontiguousarray(invL, dtype=float),
            np.ascontiguousarray(invD, dtype=float),
            np.ascontiguousarray(invZT, dtype=float),
            bool(use_ils),
            int(max_backtrack),
            out_arr,
        )
        return out_arr

    for p in range(dphase.shape[0]):
        afloat = dphase[p] / (-2.0 * np.pi)
        zfloat = Zt @ afloat
        chi2, zboot = _py_chistart(L, D, zfloat, 1)
        if use_ils:
            zcheck, _, ierr, brk = _py_lsearch(zfloat, invL, invD, chi2, max_backtrack)
            if brk:
                zcheck = zboot
        else:
            zcheck = zboot
        out_arr[p] = np.floor(invZT @ zcheck + 0.5)
    return out_arr
