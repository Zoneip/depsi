"""Per-arc estimation: integer least squares / bootstrapping + variance factors.

Port of ``ps_ils_and_bootstrap_estimation.m`` and the model-selection loop of
``ps_network.m``.  The heavy numerical part (ambiguity resolution) lives in
:mod:`depsipy.core.lambda_ils`; this module builds the geometry shared by all
arcs of one model and evaluates the a-posteriori test quantities.

The linear model is::

    y = A a + B b + e,        A = [-2*pi*I; 0],  B = [B1; I]
    Qy = diag([Qy1, diag(std_param^2)])

with ``a`` the integer ambiguities and ``b`` the physical parameters.  All
arcs share ``B1`` (built on the *mean* height-to-phase factors, as in MATLAB),
hence ``Qahat`` and the ``(B1' Qy1^-1 B1)^-1 B1' Qy1^-1`` projector are
computed once per model.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .lambda_ils import LambdaModel
from .models import NPAR_MAX, build_functional_model, model_spec

__all__ = ["Geometry", "build_geometry", "EstimationResult", "estimate_arcs"]

_DEG = np.pi / 180.0

#: above this condition number the normal-matrix inverse is treated as singular
_COND_WARN = 1e12


def _inverse(A: np.ndarray, *, what: str = "matrix") -> np.ndarray:
    """Invert a symmetric positive (semi-)definite matrix robustly.

    Tries Cholesky; on failure (or an absurd condition number) falls back to the
    pseudo-inverse, which is what a rank-deficient functional model requires
    (e.g. a constant height-to-phase factor, where topography and the master
    atmosphere term are indistinguishable).
    """
    A = np.asarray(A, dtype=float)
    cond = np.linalg.cond(A)
    if not np.isfinite(cond) or cond > _COND_WARN:
        warnings.warn(
            f"{what} is (near) singular, cond={cond:.3g}; using the "
            "pseudo-inverse. The functional model is probably rank deficient "
            "(e.g. a height-to-phase factor that is constant or nearly so).",
            RuntimeWarning,
            stacklevel=3,
        )
        return np.linalg.pinv(A)
    try:
        c = np.linalg.cholesky(A)
    except np.linalg.LinAlgError:
        warnings.warn(
            f"{what} is not positive definite; using the pseudo-inverse "
            "(the functional model is rank deficient)",
            RuntimeWarning,
            stacklevel=3,
        )
        return np.linalg.pinv(A)
    inv = np.linalg.solve(c.T, np.linalg.solve(c, np.eye(A.shape[0])))
    return 0.5 * (inv + inv.T)


@dataclass
class Geometry:
    """Matrices shared by all arcs for one functional model."""

    model: int
    Npar: int
    par_index: tuple[int, ...]
    B1: np.ndarray  # (Nifgs, Npar)
    Qy2_diag: np.ndarray  # (Npar,)
    Qy1: np.ndarray  # (Nifgs, Nifgs)
    invQy1: np.ndarray
    Qahat: np.ndarray  # (Nifgs, Nifgs)
    rhs_fixed: np.ndarray  # (Npar, Nifgs)
    b_redundancy: int
    lambda_model: LambdaModel

    @property
    def sigma0_apriori(self):
        return None


def build_geometry(
    model: int,
    Nifgs: int,
    mean_h2ph: np.ndarray,
    Btemp: np.ndarray,
    Bdop: np.ndarray | None,
    std_param: np.ndarray,
    sig2_est: np.ndarray | None,
    *,
    m2ph: float,
    sat_vel: float,
    breakpoint: int | None = None,
    breakpoint2: int | None = None,
) -> Geometry:
    """Build the least-squares geometry of one model for the network estimation.

    ``sig2_est`` is the variance-component vector from ``ps_vce.m``: entry 0 is
    the master-atmosphere variance (zero when the master is estimated) and
    entries 1..Nifgs are the per-interferogram variances.
    """
    spec = model_spec(model)
    Npar = spec.npar
    Btemp = np.asarray(Btemp, dtype=float).ravel()
    B1, Qy2_diag = build_functional_model(
        model,
        Nifgs,
        mean_h2ph,
        Btemp,
        Bdop,
        std_param,
        m2ph=m2ph,
        sat_vel=sat_vel,
        breakpoint=breakpoint,
        breakpoint2=breakpoint2,
    )

    if sig2_est is None:
        if model == 1:  # master atmosphere stochastic
            sig2 = np.concatenate(([(np.pi * 15 * _DEG) ** 2], np.full(Nifgs, (np.pi * 20 * _DEG) ** 2)))
        else:
            sig2 = np.concatenate(([0.0], np.full(Nifgs, (np.pi * 30 * _DEG) ** 2)))
    else:
        sig2 = np.asarray(sig2_est, dtype=float).ravel()
        if sig2.size != Nifgs + 1:
            raise ValueError(f"sig2_est must have {Nifgs + 1} entries, got {sig2.size}")

    # Qy1 = 2*sig2[0] * ones + diag(2*sig2[1:])
    Qy1 = np.full((Nifgs, Nifgs), 2.0 * sig2[0])
    Qy1[np.diag_indices(Nifgs)] += 2.0 * sig2[1:]

    A1 = -2 * np.pi * np.eye(Nifgs)
    A = np.vstack([A1, np.zeros((Npar, Nifgs))])
    B = np.vstack([B1, np.eye(Npar)])
    C = np.hstack([A, B])

    Qy2 = np.diag(Qy2_diag)
    Qy = np.zeros((Nifgs + Npar, Nifgs + Npar))
    Qy[:Nifgs, :Nifgs] = Qy1
    Qy[Nifgs:, Nifgs:] = Qy2

    # N = C' Qy^-1 C and its inverse.  Building inv(Qy) explicitly and then
    # inverting C'Qy^-1C (as the MATLAB code does: inv(C'*invQy*C)) squares the
    # condition number; solving keeps the error at the level of the design.
    N = C.T @ np.linalg.solve(Qy, C)
    N = 0.5 * (N + N.T)
    Qchat = _inverse(N, what="normal matrix C'Qy^-1C")

    Qahat = Qchat[:Nifgs, :Nifgs]
    Qahat = 0.5 * (Qahat + Qahat.T)

    invQy1 = _inverse(Qy1, what="Qy1")
    Nb = B1.T @ (invQy1 @ B1)
    Nb = 0.5 * (Nb + Nb.T)
    Qbhat2 = _inverse(Nb, what="normal matrix B1'Qy1^-1B1")
    rhs_fixed = Qbhat2 @ B1.T @ invQy1

    return Geometry(
        model=model,
        Npar=Npar,
        par_index=spec.par_index,
        B1=B1,
        Qy2_diag=Qy2_diag,
        Qy1=Qy1,
        invQy1=invQy1,
        Qahat=Qahat,
        rhs_fixed=rhs_fixed,
        b_redundancy=Nifgs - Npar,
        lambda_model=LambdaModel(Qahat),
    )


# ---------------------------------------------------------------------------
# per-arc evaluation
# ---------------------------------------------------------------------------


@dataclass
class EstimationResult:
    """Per-arc estimation output (arrays have one row per arc)."""

    acheck: np.ndarray  # (Narcs, Nifgs) integer ambiguities
    param: np.ndarray  # (Narcs, Npar) estimated parameters
    varfac: np.ndarray  # (Narcs,) a-posteriori variance factor
    ens_coh: np.ndarray  # (Narcs,) ensemble coherence
    model: np.ndarray  # (Narcs,) model index actually used
    althyp: np.ndarray  # (Narcs,) alternative hypothesis used (breakpoint - 1)
    bcheck: np.ndarray  # (Narcs, Npar) redundant but handy for diagnostics

    def __len__(self) -> int:
        return self.acheck.shape[0]


def evaluate(geometry: Geometry, dphase: np.ndarray, *, use_ils: bool = True) -> EstimationResult:
    """Evaluate one model on all arcs (one pass, no model selection)."""
    dphase = np.asarray(dphase, dtype=float)
    narcs, nifgs = dphase.shape
    g = geometry
    if nifgs != g.Qy1.shape[0]:
        raise ValueError(f"dphase has {nifgs} interferograms, geometry expects {g.Qy1.shape[0]}")

    acheck = g.lambda_model.resolve(dphase, use_ils=use_ils)
    y_unwrap = dphase + 2 * np.pi * acheck
    bcheck = y_unwrap @ g.rhs_fixed.T  # (Narcs, Npar)
    ycheck = bcheck @ g.B1.T
    echeck = y_unwrap - ycheck

    OMT = np.einsum("ij,jk,ik->i", echeck, g.invQy1, echeck)
    with np.errstate(divide="ignore", invalid="ignore"):
        varfac = OMT / g.b_redundancy
    ens_coh = np.abs(np.exp(1j * echeck).mean(axis=1))

    return EstimationResult(
        acheck=acheck,
        param=bcheck,
        varfac=varfac,
        ens_coh=ens_coh,
        model=np.full(narcs, g.model, dtype=float),
        althyp=np.full(narcs, np.nan),
        bcheck=bcheck,
    )


def estimate_arcs(
    dphase: np.ndarray,
    Nifgs: int,
    mean_h2ph: np.ndarray,
    Btemp: np.ndarray,
    Bdop: np.ndarray | None,
    std_param: np.ndarray,
    sig2_est: np.ndarray | None,
    models: list[int] | tuple[int, ...],
    *,
    m2ph: float,
    sat_vel: float,
    method: str = "ils",
    breakpoint: int | None = None,
    breakpoint2: int | None = None,
    varfac_switch: float = 1.0,
) -> EstimationResult:
    """Estimate all arcs, evaluating alternative models until the OMT passes.

    This reproduces the control flow of ``ps_ils_and_bootstrap_estimation.m``:
    models are tried in order; an arc is retried with the next model while its
    a-posteriori variance factor exceeds ``varfac_switch`` (1.0, i.e. the
    overall model test at 1 sigma).
    """
    if method not in ("ils", "boot"):
        raise ValueError(f"unknown ps_method {method!r} (expected 'ils' or 'boot')")
    models = [int(m) for m in np.atleast_1d(models)]
    if not models:
        raise ValueError("no model specified")

    dphase = np.asarray(dphase, dtype=float)
    narcs = dphase.shape[0]
    out = EstimationResult(
        acheck=np.full((narcs, Nifgs), np.nan),
        param=np.full((narcs, NPAR_MAX), np.nan),
        varfac=np.full(narcs, np.inf),
        ens_coh=np.full(narcs, np.nan),
        model=np.full(narcs, np.nan),
        althyp=np.full(narcs, np.nan),
        bcheck=np.full((narcs, NPAR_MAX), np.nan),
    )

    todo = np.arange(narcs)
    for model in models:
        if todo.size == 0:
            break
        geom = build_geometry(
            model,
            Nifgs,
            mean_h2ph,
            Btemp,
            Bdop,
            std_param,
            sig2_est,
            m2ph=m2ph,
            sat_vel=sat_vel,
            breakpoint=breakpoint,
            breakpoint2=breakpoint2,
        )
        res = evaluate(geom, dphase[todo], use_ils=(method == "ils"))
        improved = res.varfac < out.varfac[todo]
        idx = todo[improved]
        if idx.size:
            out.acheck[idx] = res.acheck[improved]
            out.param[np.ix_(idx, geom.par_index - 1)] = res.param[improved]
            out.varfac[idx] = res.varfac[improved]
            out.ens_coh[idx] = res.ens_coh[improved]
            out.model[idx] = model
            out.althyp[idx] = (breakpoint - 1) if breakpoint is not None else np.nan
            out.bcheck[np.ix_(idx, geom.par_index - 1)] = res.bcheck[improved]
        # the original stops testing once the OMT passes (varfac <= 1)
        todo = todo[~(res.varfac <= varfac_switch)]
    return out
