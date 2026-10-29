"""Functional models for PS estimation (port of ``ps_model_definitions.m``).

A model describes the *phase* as a linear function of ``Npar`` parameters::

    phase(Nifgs) = B1 @ param

with the first parameter always the topographic height.  Additionally each
model carries a *pseudo-observation* per parameter (the ``std_param`` entries)
which regularises the estimation; these enter the stochastic model as the
diagonal block ``Qy2``.

The MATLAB original stores the parameter list in a global and mutates
``althyp_index`` inside the function, which makes the behaviour depend on the
number of previously evaluated models.  Here the breakpoints are explicit
arguments, which is what the original *intends* (see the two-loop control flow
in ``ps_ils_and_bootstrap_estimation.m``); see ``docs/parity.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "NPAR_MAX",
    "FIRST_DEFO_INDEX",
    "PAR_NAMES",
    "PAR_UNITS",
    "PAR_CAXIS",
    "ModelSpec",
    "model_spec",
    "npar",
    "build_functional_model",
    "covariance_index",
    "defo_index",
    "annotation",
]

#: maximum number of parameters of any model (``Npar_max`` in the MATLAB code)
NPAR_MAX = 10
#: indices >= this are deformation parameters
FIRST_DEFO_INDEX = 4

#: 1-based parameter identifiers, exactly as in ``ps_model_definitions.m``
PAR_NAMES = {
    1: "topo",
    2: "master_noise",
    3: "subpixel",
    4: "linear1",
    5: "linear2",
    6: "linear3",
    7: "quadratic",
    8: "cubic",
    9: "periodic1",
    10: "periodic2",
}
PAR_UNITS = {
    1: "m",
    2: "m",
    3: "m",
    4: "m/y",
    5: "m/y",
    6: "m/y",
    7: "m/y^2",
    8: "m/y^3",
    9: "-",
    10: "-",
}
PAR_CAXIS = {
    1: (-20.0, 20.0),
    2: (-0.01, 0.01),
    3: (-5.0, 5.0),
    4: (-0.01, 0.01),
    5: (-0.01, 0.01),
    6: (-0.01, 0.01),
    7: (-0.01, 0.01),
    8: (-0.01, 0.01),
    9: (-0.05, 0.05),
    10: (-0.05, 0.05),
}


@dataclass(frozen=True)
class ModelSpec:
    """A named functional model: the 1-based parameter indices it contains."""

    index: int
    par_index: tuple[int, ...]
    description: str

    @property
    def npar(self) -> int:
        return len(self.par_index)

    @property
    def defo_index(self) -> tuple[int, ...]:
        """0-based positions (within the parameter vector) of deformation parameters."""
        return tuple(i for i, p in enumerate(self.par_index) if p >= FIRST_DEFO_INDEX)


_MODELS: dict[int, tuple[tuple[int, ...], str]] = {
    1: ((1, 4), "linear for periodogram"),
    2: ((1, 2, 4), "linear"),
    3: ((1, 2, 4, 7), "quadratic"),
    4: ((1, 2, 4, 7, 8), "cubic"),
    5: ((1, 2, 4, 9, 10), "linear+periodic"),
    6: ((1, 2, 4, 7, 9, 10), "quadratic+periodic"),
    7: ((1, 2, 4, 7, 8, 9, 10), "cubic+periodic"),
    8: ((1, 2, 3, 4), "linear+subpixel"),
    9: ((1, 2, 3, 4, 7), "quadratic+subpixel"),
    10: ((1, 2, 3, 4, 7, 8), "cubic+subpixel"),
    11: ((1, 2, 3, 4, 9, 10), "linear+periodic+subpixel"),
    12: ((1, 2, 3, 4, 7, 9, 10), "quadratic+periodic+subpixel"),
    13: ((1, 2, 3, 4, 7, 8, 9, 10), "cubic+periodic+subpixel"),
    14: ((1, 2, 4, 5), "double linear"),
    15: ((1, 2, 4, 5, 7), "double quadratic"),
    16: ((1, 2, 4, 5, 7, 8), "double cubic"),
    17: ((1, 2, 4, 5, 9, 10), "double linear+periodic"),
    18: ((1, 2, 4, 5, 7, 9, 10), "double quadratic+periodic"),
    19: ((1, 2, 4, 5, 7, 8, 9, 10), "double cubic+periodic"),
    20: ((1, 2, 3, 4, 5), "double linear+subpixel"),
    21: ((1, 2, 3, 4, 5, 7), "double quadratic+subpixel"),
    22: ((1, 2, 3, 4, 5, 7, 8), "double cubic+subpixel"),
    23: ((1, 2, 3, 4, 5, 9, 10), "double linear+periodic+subpixel"),
    24: ((1, 2, 3, 4, 5, 7, 9, 10), "double quadratic+periodic+subpixel"),
    25: ((1, 2, 3, 4, 5, 7, 8, 9, 10), "double cubic+periodic+subpixel"),
    26: ((1, 2, 4, 5, 6), "triple"),
    27: ((1, 2, 3, 4, 5, 6), "triple+subpixel"),
}

#: models that use a single breakpoint (14..25) / two breakpoints (26..27)
SINGLE_BREAKPOINT_MODELS = range(14, 26)
DOUBLE_BREAKPOINT_MODELS = range(26, 28)


def model_spec(model: int) -> ModelSpec:
    """Return the :class:`ModelSpec` of model ``model`` (1..27)."""
    idx = int(np.asarray(model).ravel()[0])
    if idx not in _MODELS:
        raise ValueError(
            f"model {idx} is not implemented (valid: 1..27, see ps_model_definitions.m)"
        )
    par, desc = _MODELS[idx]
    return ModelSpec(index=idx, par_index=par, description=desc)


def npar(model: int) -> int:
    return model_spec(model).npar


def covariance_index(spec: ModelSpec) -> np.ndarray:
    """Reproduce ``covar_index`` from ``ps_model_definitions.m`` (1-based values)."""
    par = np.asarray(spec.par_index, dtype=int)
    np_ = par.size
    out: list[int] = []
    for v in range(np_):
        offs = np.arange(NPAR_MAX - par[v] + 1, NPAR_MAX)  # sum(...) in MATLAB
        out.extend((par[v:] + int(offs.sum())).tolist())
    return np.asarray(out, dtype=int)


def defo_index(model: int) -> tuple[int, ...]:
    return model_spec(model).defo_index


def annotation(model: int) -> tuple[list[str], list[str]]:
    """Names and units of the parameters of ``model`` (length ``NPAR_MAX``)."""
    spec = model_spec(model)
    names = ["not_used"] * NPAR_MAX
    units = ["not_used"] * NPAR_MAX
    for k, par in enumerate(spec.par_index):
        names[k] = PAR_NAMES[par]
        units[k] = PAR_UNITS[par]
    return names, units


# ---------------------------------------------------------------------------
# design matrix construction
# ---------------------------------------------------------------------------


def _linear_basis(Btemp: np.ndarray, bp1: int | None, bp2: int | None) -> np.ndarray:
    """First deformation basis function (``linear1``), 0-based ``bp`` indices."""
    n = Btemp.size
    if bp1 is None:
        return Btemp.copy()
    b1 = Btemp[bp1]
    if bp2 is None:
        if b1 > 0:
            out = Btemp.copy()
            out[bp1:] = b1
            return out
        if b1 < 0:
            out = np.zeros(n)
            out[:bp1] = Btemp[:bp1] - b1
            return out
        raise ValueError("breakpoint has a zero temporal baseline")
    b2 = Btemp[bp2]
    if b1 >= 0 and b2 >= 0:
        out = Btemp.copy()
        out[bp1:] = b1
        return out
    if b1 < 0:
        out = np.zeros(n)
        out[:bp1] = Btemp[:bp1] - b1
        return out
    raise ValueError("bp2 must be larger than bp1")


def _quadratic_basis(Btemp: np.ndarray, bp1: int | None, bp2: int | None) -> np.ndarray:
    """Second deformation basis function (``linear2``)."""
    n = Btemp.size
    if bp1 is None:
        raise ValueError("a deformation model with a second linear term needs a breakpoint")
    b1 = Btemp[bp1]
    if bp2 is None:
        if b1 > 0:
            out = np.zeros(n)
            out[bp1:] = Btemp[bp1:] - b1
            return out
        if b1 < 0:
            out = Btemp.copy()
            out[:bp1] = b1
            return out
        raise ValueError("breakpoint has a zero temporal baseline")
    b2 = Btemp[bp2]
    if b1 >= 0 and b2 >= 0:
        out = np.zeros(n)
        out[bp1:bp2] = Btemp[bp1:bp2] - b1
        out[bp2:] = b2 - b1
        return out
    if b1 < 0 and b2 >= 0:
        out = Btemp.copy()
        out[:bp1] = b1
        out[bp1:bp2] = Btemp[bp1:bp2]
        out[bp2:] = b2
        return out
    if b1 < 0 and b2 < 0:
        out = Btemp.copy()
        out[:bp1] = b1 - b2
        out[bp1:bp2] = Btemp[bp1:bp2] - b2
        out[bp2:] = 0.0
        return out
    raise ValueError("bp2 must be larger than bp1")


def _cubic_basis(Btemp: np.ndarray, bp1: int | None, bp2: int | None) -> np.ndarray:
    """Third deformation basis function (``linear3``)."""
    n = Btemp.size
    if bp1 is None or bp2 is None:
        raise ValueError("a deformation model with a third linear term needs two breakpoints")
    b1, b2 = Btemp[bp1], Btemp[bp2]
    out = np.zeros(n)
    if b1 >= 0 and b2 >= 0:
        out[bp2:] = Btemp[bp2:] - b2
        return out
    if b1 < 0 and b2 >= 0:
        out[bp2:] = Btemp[bp2:] - b2
        return out
    if b1 < 0 and b2 < 0:
        out[:bp2] = b2
        out[bp2:] = Btemp[bp2:]
        return out
    raise ValueError("bp2 must be larger than bp1")


def build_functional_model(
    model: int,
    Nifgs: int,
    h2ph: np.ndarray,
    Btemp: np.ndarray,
    Bdop: np.ndarray | None,
    std_param: np.ndarray,
    *,
    m2ph: float,
    sat_vel: float,
    breakpoint: int | None = None,
    breakpoint2: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the design matrix ``B1`` and the pseudo-observation variances.

    Parameters
    ----------
    model
        Model index 1..27 (see :func:`model_spec`).
    Nifgs
        Number of interferograms.
    h2ph
        Height-to-phase factors, shape ``(Nifgs,)`` or ``(Narcs, Nifgs)``; the
        column mean is used, exactly as in the MATLAB code (``mean_h2ph``).
    Btemp
        Temporal baselines in years, shape ``(Nifgs,)``.
    Bdop
        Doppler baselines in Hz (only used by models with a subpixel term).
    std_param
        1..7 standard deviations of the pseudo-observations.
    m2ph
        metres-to-phase constant ``-4*pi/lambda``.
    sat_vel
        Satellite velocity, used by the subpixel term.
    breakpoint, breakpoint2
        1-based interferogram indices (MATLAB convention: the break lies
        between ``breakpoint-1`` and ``breakpoint``).

    Returns
    -------
    B1 : ndarray, shape (Nifgs, Npar)
        Phase design matrix (columns in metres-equivalent phase units).
    Qy2_diag : ndarray, shape (Npar,)
        Variances of the pseudo-observations (diagonal of ``Qy2``).
    """
    spec = model_spec(model)
    Btemp = np.asarray(Btemp, dtype=float).ravel()
    if Btemp.size != Nifgs:
        raise ValueError(f"Btemp has {Btemp.size} entries, expected {Nifgs}")
    h2ph = np.asarray(h2ph, dtype=float)
    mean_h2ph = h2ph.mean(axis=0).ravel() if h2ph.ndim == 2 else h2ph.ravel()
    if mean_h2ph.size != Nifgs:
        raise ValueError(f"h2ph gives {mean_h2ph.size} columns, expected {Nifgs}")
    std = np.asarray(std_param, dtype=float).ravel()

    if model in SINGLE_BREAKPOINT_MODELS and breakpoint is None:
        raise ValueError(f"model {model} requires a breakpoint")
    if model in DOUBLE_BREAKPOINT_MODELS and (breakpoint is None or breakpoint2 is None):
        raise ValueError(f"model {model} requires two breakpoints")

    bp1 = None if breakpoint is None else int(breakpoint) - 1
    bp2 = None if breakpoint2 is None else int(breakpoint2) - 1
    for name, bp in (("breakpoint", bp1), ("breakpoint2", bp2)):
        if bp is not None and not (0 <= bp < Nifgs):
            raise ValueError(f"{name}={bp + 1} outside 1..{Nifgs}")

    B1 = np.zeros((Nifgs, spec.npar))
    Qy2 = np.zeros(spec.npar)

    def _std(k: int) -> float:
        if k - 1 >= std.size:
            raise ValueError(f"std_param has {std.size} entries, model uses entry {k}")
        return float(std[k - 1]) ** 2

    for col, par in enumerate(spec.par_index):
        if par == 1:
            B1[:, col] = m2ph * mean_h2ph
            Qy2[col] = _std(1)
        elif par == 2:
            B1[:, col] = m2ph
            Qy2[col] = _std(2)
        elif par == 3:
            if Bdop is None:
                raise ValueError("model requires Doppler baselines (Bdop)")
            B1[:, col] = (2 * np.pi / sat_vel) * np.asarray(Bdop, dtype=float).ravel()
            Qy2[col] = _std(3)
        elif par == 4:
            B1[:, col] = m2ph * _linear_basis(Btemp, bp1, bp2)
            Qy2[col] = _std(4)
        elif par == 5:
            B1[:, col] = m2ph * _quadratic_basis(Btemp, bp1, bp2)
            Qy2[col] = _std(4)
        elif par == 6:
            B1[:, col] = m2ph * _cubic_basis(Btemp, bp1, bp2)
            Qy2[col] = _std(4)
        elif par == 7:
            B1[:, col] = m2ph * Btemp**2
            Qy2[col] = _std(5)
        elif par == 8:
            B1[:, col] = m2ph * Btemp**3
            Qy2[col] = _std(6)
        elif par == 9:
            B1[:, col] = m2ph * np.sin(2 * np.pi * Btemp)
            Qy2[col] = _std(7)
        elif par == 10:
            B1[:, col] = m2ph * (np.cos(2 * np.pi * Btemp) - 1.0)
            Qy2[col] = _std(7)
        else:  # pragma: no cover - guarded by model_spec
            raise ValueError(f"unknown parameter {par}")

    return B1, Qy2
