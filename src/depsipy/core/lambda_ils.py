"""LAMBDA integer least squares + bootstrapping for PS ambiguity resolution.

This is a port of the Delft LAMBDA routines used by DePSI
(``decorrel_freek.m``, ``ldldecom.m``, ``ps_chistart.m``, ``bk_lsearch_count.m``)
together with the per-arc driver of ``ps_ils_and_bootstrap_estimation.m``.

Design notes
------------
The original MATLAB code spends most of its time in a ``parfor`` over arcs with
a *scalar* loop of the integer search inside.  Here the same algorithm is
implemented as a compiled (numba) kernel over all arcs at once, with the initial
search-ellipsoid size reused across arcs; the expensive part of the estimator --
the conditional rounding for the bootstrapped candidate -- is exactly the same
arithmetic, so results agree with MATLAB.

If numba is unavailable the module falls back to a vectorised numpy
implementation which is ~20x slower but numerically identical.
"""

from __future__ import annotations

import numpy as np

from . import _kernels as K

__all__ = [
    "ldldecom",
    "decorrel",
    "chistart",
    "lsearch",
    "LambdaModel",
    "resolve_ambiguities",
]


from ._reference import chistart, decorrel, ldldecom, lsearch  # noqa: F401  (re-exported)

# ---------------------------------------------------------------------------
# the reusable per-model LAMBDA state
# ---------------------------------------------------------------------------


class LambdaModel:
    """Pre-computed LAMBDA quantities for one functional/stochastic model.

    All arcs share ``Qahat`` (it only depends on the design matrix and the
    variance components), so the decorrelation is done once.
    """

    __slots__ = ("Z", "L", "D", "invZ", "invL", "invD", "invLtDL", "n", "max_backtrack")

    def __init__(self, Qahat: np.ndarray, *, max_backtrack: int | None = None, compat: bool = False):
        Qahat = np.asarray(Qahat, dtype=float)
        self.n = Qahat.shape[0]
        # MATLAB uses n**3 as a hard guard and falls back to the bootstrapped
        # candidate when the search exceeds it, which is *not* optimal; the
        # default here gives the branch-and-bound room to finish (see docs).
        if max_backtrack is None:
            max_backtrack = self.n**3 if compat else max(100_000, (self.n**3) * 1000)
        self.max_backtrack = int(max_backtrack)
        self.Z, self.L, self.D = decorrel(Qahat)
        self.invZ = np.linalg.inv(self.Z.astype(float))
        self.invL = np.linalg.inv(self.L)
        self.invD = 1.0 / self.D
        self.invLtDL = np.linalg.inv(self.L.T @ np.diag(self.D) @ self.L)

    # -- single arc --------------------------------------------------------
    def resolve_one(self, dphase_arc: np.ndarray, use_ils: bool = True):
        """Resolve one arc.  Returns ``(acheck, zcheck, boot, chi2, break_flag)``."""
        afloat = np.asarray(dphase_arc, dtype=float) / (-2 * np.pi)
        zfloat = self.Z.T @ afloat
        chi2, zboot = chistart(self.L, self.D, zfloat, 1)
        if use_ils:
            zcheck, _, ierr, brk = lsearch(zfloat, self.invL, self.invD, chi2, self.max_backtrack)
            if brk:
                zcheck = zboot
        else:
            zcheck, brk = zboot, False
        acheck = np.round(self.invZ.T @ zcheck)
        return acheck, zcheck, zboot, chi2, brk

    # -- many arcs ---------------------------------------------------------
    def resolve(self, dphase: np.ndarray, use_ils: bool = True) -> np.ndarray:
        """Resolve all arcs; ``dphase`` is ``(Narcs, Nifgs)``. Returns ``acheck``."""
        dphase = np.ascontiguousarray(np.asarray(dphase, dtype=float))
        if dphase.ndim != 2 or dphase.shape[1] != self.n:
            raise ValueError(f"dphase must be (Narcs, {self.n}), got {dphase.shape}")
        out = np.full(dphase.shape, np.nan)
        K.resolve_all(
            dphase,
            self.Z.T.astype(float),
            self.invLtDL,
            self.L,
            self.D,
            self.invL,
            self.invD,
            self.invZ.T.astype(float),
            bool(use_ils),
            int(self.max_backtrack),
            out,
        )
        return out


def resolve_ambiguities(Qahat, dphase, use_ils=True) -> np.ndarray:
    """Convenience wrapper: decorrelate ``Qahat`` and resolve all arcs."""
    return LambdaModel(Qahat).resolve(dphase, use_ils=use_ils)
