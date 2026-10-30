"""End-to-end processing chain (groups 2..8 of ``depsi.m``).

The MATLAB driver keeps its state in ``<project_id>_project.mat`` and reloads it
between processing groups.  Here the same stages are plain functions over a
:class:`Project` context, so a run is one Python process and every intermediate
is inspectable:

====================  ==============================================  ==========
group                 stage                                           module
====================  ==============================================  ==========
2                     PSC / PSP selection                             selection
3                     VCE + network + estimation + unwrapping         network / core / unwrapping
3a                    trend removal (plane per IFG, with fallback)    this module
4                     second network on detrended data, atmosphere    network / filtering
5/6                   densification of PSP                            densification
7                     spatio-temporal consistency + results           stc / output
8                     temporal noise filtering                        noise_filter
====================  ==============================================  ==========

Files are written with the DePSI names/layouts (``_psc_2orig_sel*.raw``,
``_psc_results_sel*.raw``, ``_psc_atmo_filt_sel*.raw`` ...), so the MATLAB side
and the readiness tools can read them.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .core.estimation import covariance_qy1, estimate_arcs, stochastic_sig2
from .core.models import NPAR_MAX, build_functional_model, model_spec
from .io.params import Params, load_params
from .io.project import InputStack, read_input_file, sensor_constants
from .io.raw import write_raw
from .steps.densification import DensificationResult, densify_points
from .io.doris import read_doris_metadata
from .steps.filtering import separate_atmosphere
from .steps.geocoding import lph2xyz, orbit_state, xyz_to_latlonh
from .steps.network import ArcNetwork, estimate_variance_components, form_network
from .steps.noise_filter import filter_deformation
from .steps.output import write_point_csv, write_point_stack, write_resfile, write_summary
from .steps.selection import SelectionResult, select_psc_psp, write_psc_files
from .steps.stc import spatiotemporal_consistency
from .steps.unwrapping import UnwrapResult, unwrap_network

__all__ = [
    "Project",
    "NetworkResult",
    "PointFit",
    "DetrendResult",
    "load_project",
    "run_selection",
    "run_network",
    "fit_points",
    "run_detrend",
    "run_densification",
    "run_filtering",
    "run_stc",
    "run_noise_filter",
    "run_output",
    "run_pipeline",
]


# ---------------------------------------------------------------------------
# project context
# ---------------------------------------------------------------------------


@dataclass
class Project:
    """Everything a stage needs: parameters, stack geometry and file lists."""

    params: Params
    stack: InputStack
    workdir: Path
    project_id: str
    crop: tuple[int, int, int, int]
    Nlines: int
    Npixels: int
    az_spacing: float
    r_spacing: float
    lam: float
    vel: float
    m2ph: float
    slc_paths: list[str]
    ifg_paths: list[str]
    h2ph_paths: list[str]
    slc_selection: np.ndarray
    master_res: str | None = None

    @property
    def Nifgs(self) -> int:
        return self.stack.Nifgs

    @property
    def Nslc(self) -> int:
        return self.stack.Nslc

    def path(self, suffix: str) -> Path:
        """Path of a DePSI product ``<project_id><suffix>`` in the work directory."""
        return self.workdir / f"{self.project_id}{suffix}"


def load_project(param_file, *, workdir=None) -> Project:
    """Read the parameter file and the input file, resolve every raster path."""
    params = load_params(param_file)
    param_dir = Path(param_file).resolve().parent
    stack = read_input_file(params.str("input_file"))
    if workdir is None:
        workdir = param_dir
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    crop_final = stack.crop_final
    if params.has("crop") and not params.is_empty("crop"):
        requested = tuple(int(v) for v in np.asarray(params.get("crop")).ravel()[:4])
        crop = (
            max(requested[0], crop_final[0]),
            min(requested[1], crop_final[1]),
            max(requested[2], crop_final[2]),
            min(requested[3], crop_final[3]),
        )
        if crop[1] < crop[0] or crop[3] < crop[2]:
            raise ValueError(f"crop {requested} does not intersect the input crop {crop_final}")
    else:
        crop = crop_final

    lam, vel, m2ph = sensor_constants(params.str("sensor"))
    slc_paths = stack.slc_paths()
    ifg_paths = [stack.resolve(f) for f in stack.ifg]
    h2ph_paths = [stack.resolve(f) for f in stack.h2ph]

    if params.has("slc_selection_input") and not params.is_empty("slc_selection_input"):
        slc_selection = params.ints("slc_selection_input")
    else:
        slc_selection = np.arange(1, stack.Nslc + 1)

    master_res = None
    if params.has("master_res") and not params.is_empty("master_res"):
        master_res = params.str("master_res")
        candidate = Path(master_res)
        master_res = str(candidate if candidate.is_absolute() else param_dir / candidate)

    return Project(
        params=params,
        stack=stack,
        workdir=workdir,
        project_id=params.str("project_id"),
        crop=crop,
        Nlines=crop[1] - crop[0] + 1,
        Npixels=crop[3] - crop[2] + 1,
        az_spacing=params.float("az_spacing"),
        r_spacing=params.float("r_spacing"),
        lam=lam,
        vel=vel,
        m2ph=m2ph,
        slc_paths=slc_paths,
        ifg_paths=ifg_paths,
        h2ph_paths=h2ph_paths,
        slc_selection=slc_selection,
        master_res=master_res,
    )


# ---------------------------------------------------------------------------
# group 2: selection
# ---------------------------------------------------------------------------


def run_selection(project: Project, *, block_rows: int | None = None) -> SelectionResult:
    """Select PSC/PSP and write ``_psc_2orig_sel*`` / ``_psp_2orig_sel*``."""
    p = project.params
    ref_cn = None if p.is_empty("ref_cn") else p.floats("ref_cn")
    result = select_psc_psp(
        slc_files=project.slc_paths,
        ifg_files=project.ifg_paths,
        h2ph_files=project.h2ph_paths,
        Nlines=project.Nlines,
        Npixels=project.Npixels,
        az_spacing=project.az_spacing,
        r_spacing=project.r_spacing,
        crop=project.crop,
        crop_in=project.stack.crop_in,
        slc_selection=project.slc_selection,
        calfactors=np.ones(project.slc_selection.size),
        gridsize=p.float("psc_selection_gridsize"),
        psc_threshold=p.float("psc_threshold"),
        psp_threshold1=p.float("psp_threshold1"),
        psp_threshold2=None if p.is_empty("psp_threshold2") else p.float("psp_threshold2"),
        psc_distribution=p.str("psc_distribution", "nonuniform"),
        do_apriori_sidelobe_mask=p.str("do_apriori_sidelobe_mask", "yes"),
        livetime_threshold=p.float("livetime_threshold", 0.2),
        peak_tolerance=p.float("peak_tolerance", 0.99),
        Npsc_selections=p.int("Npsc_selections", 1),
        ref_cn=ref_cn,
        workspace=project.workdir,
        block_rows=block_rows,
    )
    write_psc_files(result, project.workdir, project.project_id, project.Nifgs)
    if result.mean_amp is not None:
        write_raw(project.path("_mrm.raw"), result.mean_amp, "float32")
    if result.amp_disp is not None:
        write_raw(project.path("_amp_disp.raw"), result.amp_disp, "float32")
    if result.side_lobe is not None:
        write_raw(project.path("_side_lobe_mask.raw"), result.side_lobe.astype(np.float32), "float32")
    return result


# ---------------------------------------------------------------------------
# group 3/4: network estimation and unwrapping
# ---------------------------------------------------------------------------


@dataclass
class PointFit:
    """Per-PS fit of the final model to the unwrapped phase."""

    param: np.ndarray  # (Nps, NPAR_MAX)
    ens_coh: np.ndarray
    sig2hat: np.ndarray
    phase_res: np.ndarray  # (Nps, Nifgs)
    defo: np.ndarray  # (Nps, Nifgs) metres, model deformation + residual
    model: int


@dataclass
class NetworkResult:
    """Output of one ``ps_network``-equivalent pass (VCE -> net -> est -> unwrap)."""

    z: int
    arcs: ArcNetwork
    estimate: object  # core.estimation.EstimationResult
    unwrap: UnwrapResult
    sig2_est: np.ndarray
    model: int
    final_model: int
    psc: PointFit
    Nnodes: int

    @property
    def acheck(self) -> np.ndarray:
        return self.unwrap.acheck


def fit_points(
    phase: np.ndarray,
    h2ph: np.ndarray,
    model: int,
    Btemp: np.ndarray,
    Bdop: np.ndarray | None,
    std_param: np.ndarray,
    sig2_est: np.ndarray | None,
    *,
    m2ph: float,
    sat_vel: float,
    breakpoint: int | None = None,
    breakpoint2: int | None = None,
) -> PointFit:
    """Least-squares fit of ``model`` to every PS time series (``ps_spatial_unwrapping``).

    The design matrix has the height-to-phase column ``m2ph*h2ph_i`` (one per PS)
    and otherwise fixed columns, so the normal equations are built for all points
    at once.  A point whose normal matrix is singular falls back to the
    pseudo-inverse, exactly like the per-arc geometry does.
    """
    phase = np.atleast_2d(np.asarray(phase, dtype=float))
    h2ph = np.atleast_2d(np.asarray(h2ph, dtype=float))
    if phase.shape != h2ph.shape:
        raise ValueError("phase and h2ph must have the same shape")
    n_ps, n_ifg = phase.shape
    spec = model_spec(model)
    npar = spec.npar

    # Only ``Qy1`` and the fixed (h2ph-free) columns are shared; the topo column
    # is replaced per point.  A non-constant dummy h2ph is used for the design so
    # a constant one cannot make the topo column identical to the master-noise
    # column (model 2).
    dummy_h2ph = np.linspace(0.5, 1.5, n_ifg)
    B0, _ = build_functional_model(
        model, n_ifg, dummy_h2ph, Btemp, Bdop, std_param,
        m2ph=m2ph, sat_vel=sat_vel, breakpoint=breakpoint, breakpoint2=breakpoint2,
    )
    topo_col = [i for i, par in enumerate(spec.par_index) if par == 1]
    B = np.broadcast_to(B0, (n_ps, n_ifg, npar)).copy()
    if topo_col:
        B[:, :, topo_col[0]] = m2ph * h2ph
    invQy = np.linalg.inv(covariance_qy1(stochastic_sig2(model, n_ifg, sig2_est)))

    normal = np.einsum("pni,nm,pmj->pij", B, invQy, B)
    normal = 0.5 * (normal + np.swapaxes(normal, -1, -2))
    rhs = np.einsum("pni,nm,pm->pi", B, invQy, phase)

    x = np.full((n_ps, npar), np.nan)
    singular = np.zeros(n_ps, dtype=bool)
    try:
        x[:] = np.linalg.solve(normal, rhs[..., None])[..., 0]
    except np.linalg.LinAlgError:
        for k in range(n_ps):
            try:
                x[k] = np.linalg.solve(normal[k], rhs[k])
            except np.linalg.LinAlgError:
                x[k] = np.linalg.pinv(normal[k]) @ rhs[k]
                singular[k] = True

    yhat = np.einsum("pni,pi->pn", B, x)
    e = phase - yhat
    omt = np.einsum("pn,nm,pm->p", e, invQy, e)
    redun = max(n_ifg - npar, 1)
    sig2hat = omt / redun
    ens_coh = np.abs(np.exp(1j * e).mean(axis=1))

    defo = yhat * 0.0
    if spec.defo_index:
        defo = np.einsum("pni,pi->pn", B[:, :, spec.defo_index], x[:, spec.defo_index])
    defo = (defo + e) / m2ph

    param = np.full((n_ps, NPAR_MAX), np.nan)
    param[:, np.asarray(spec.par_index) - 1] = x

    return PointFit(param=param, ens_coh=ens_coh, sig2hat=sig2hat, phase_res=e, defo=defo, model=int(model))


def run_network(
    project: Project,
    selection: SelectionResult,
    *,
    z: int = 0,
    phase: np.ndarray | None = None,
    h2ph: np.ndarray | None = None,
    model: int | None = None,
    final_model: int | None = None,
    ref_cn=None,
) -> NetworkResult:
    """One network pass: VCE, arc estimation, spatial unwrapping, per-PS fit."""
    p = project.params
    cand = selection.psc[z]
    if not cand:
        raise ValueError(f"selection {z + 1} contains no PSC")
    az = np.array([c.az for c in cand], dtype=float)
    r = np.array([c.r for c in cand], dtype=float)
    if phase is None:
        phase = np.stack([c.phase for c in cand])
    if h2ph is None:
        h2ph = np.stack([c.h2ph for c in cand])
    phase = np.asarray(phase, dtype=float)
    h2ph = np.asarray(h2ph, dtype=float)

    psc_models = p.ints("psc_model", 1).tolist()
    final = int(p.ints("final_model", 1)[0]) if final_model is None else int(final_model)
    psc_model = int(model if model is not None else (psc_models[0] if psc_models else 1))
    method = p.str("ps_method", "ils")
    std_param = p.floats("std_param")
    breakpoint = None if p.is_empty("breakpoint") else p.int("breakpoint")
    breakpoint2 = None if p.is_empty("breakpoint2") else p.int("breakpoint2")

    net = form_network(
        az, r, phase, h2ph,
        az_spacing=project.az_spacing, r_spacing=project.r_spacing,
        max_arc_length=p.float("max_arc_length", 1000.0),
        method=p.str("network_method", "delaunay"),
        Ncon=p.int("Ncon", 8),
    )
    sig2 = estimate_variance_components(
        net.arcs, net.dphase, net.dh2ph,
        Nifgs=project.Nifgs, Btemp=project.stack.Btemp, Bdop=project.stack.Bdop,
        std_param=std_param, model=psc_model, m2ph=project.m2ph, sat_vel=project.vel,
        breakpoint=breakpoint, breakpoint2=breakpoint2, method=method,
    )
    est = estimate_arcs(
        net.dphase, project.Nifgs, net.dh2ph, project.stack.Btemp, project.stack.Bdop,
        std_param, sig2, psc_models, m2ph=project.m2ph, sat_vel=project.vel,
        method=method, breakpoint=breakpoint, breakpoint2=breakpoint2,
        varfac_switch=p.float("varfac_switch", 1.0),
    )
    unw = unwrap_network(
        net.arcs, net.dphase, est.acheck, est.varfac, az.size,
        ref_cn=ref_cn,
        weighted=p.str("weighted_unwrap", "yes").lower() == "yes",
        varfac_threshold=p.float("varfac_threshold", 3.0),
        ens_coh=est.ens_coh,
        ens_coh_threshold=p.float("ens_coh_threshold", 0.9),
        method=method,
        Nlines=project.Nlines, Npixels=project.Npixels,
        az_spacing=project.az_spacing, r_spacing=project.r_spacing,
        ps_az=az, ps_r=r,
    )

    # per-PS final model fit on the unwrapped phase (finite points only)
    psc_fit = fit_points(
        unw.phase_unw, h2ph, final,
        project.stack.Btemp, project.stack.Bdop, std_param, sig2,
        m2ph=project.m2ph, sat_vel=project.vel,
        breakpoint=breakpoint, breakpoint2=breakpoint2,
    )
    return NetworkResult(
        z=z, arcs=net, estimate=est, unwrap=unw, sig2_est=sig2,
        model=psc_model, final_model=final, psc=psc_fit, Nnodes=az.size,
    )


# ---------------------------------------------------------------------------
# group 3a: trend removal
# ---------------------------------------------------------------------------


@dataclass
class DetrendResult:
    param: np.ndarray  # (3, Nifgs) plane coefficients [az, r, 1]
    psc_phase: np.ndarray  # (Npsc, Nifgs) detrended, wrapped
    psp_phase: np.ndarray  # (Npsp, Nifgs) detrended, wrapped
    psc_trend: np.ndarray
    psp_trend: np.ndarray


def _wrap(phase: np.ndarray) -> np.ndarray:
    return np.mod(phase + np.pi, 2 * np.pi) - np.pi


def run_detrend(
    project: Project,
    selection: SelectionResult,
    net: NetworkResult,
    *,
    z: int = 0,
    max_residual: float = 100.0 / 180.0 * np.pi,
    reject_fraction: float = 0.5,
    max_iter: int = 500,
    min_points: int = 10,
    min_fraction: float = 0.5,
) -> DetrendResult:
    """Remove a plane per interferogram, estimated on the PSC (``ps_detrend.m``).

    The reference implementation rejects half of the residuals above 100 deg per
    iteration and stops when fewer than three points would remain, which on rough
    topography can end in an almost degenerate three-point plane with a slope two
    orders of magnitude too large.  Here the rejection stops as soon as fewer than
    ``min_fraction`` of the points (and never fewer than ``min_points``) would
    remain, keeping the last well-constrained fit; see ``docs/parity.md``.
    """
    cand = selection.psc[z]
    psc_az = np.array([c.az for c in cand], dtype=float)
    psc_r = np.array([c.r for c in cand], dtype=float)
    phase_wrapped = np.stack([c.phase for c in cand])
    phase_unw = 2 * np.pi * net.unwrap.acheck + phase_wrapped
    n_ifg = phase_unw.shape[1]
    coef = np.full((3, n_ifg), np.nan)
    A_all = np.column_stack([psc_az, psc_r, np.ones(psc_az.size)])

    for v in range(n_ifg):
        y = phase_unw[:, v]
        ok = np.isfinite(y)
        A = A_all[ok]
        yy = y[ok]
        n_start = A.shape[0]
        if n_start < 3:
            continue
        keep_min = max(int(min_points), int(min_fraction * n_start))
        xhat = np.zeros(3)
        for _ in range(int(max_iter)):
            xhat = np.linalg.lstsq(A, yy, rcond=None)[0]
            ehat = yy - A @ xhat
            if np.max(np.abs(ehat)) <= max_residual:
                break
            n_rej = max(1, int(np.ceil(reject_fraction * np.sum(np.abs(ehat) > max_residual))))
            if A.shape[0] - n_rej < keep_min:
                break
            drop = np.argsort(np.abs(ehat))[-n_rej:]
            A = np.delete(A, drop, axis=0)
            yy = np.delete(yy, drop)
        coef[:, v] = xhat
    if np.any(~np.isfinite(coef)):
        raise ValueError("trend estimation failed for some interferograms")

    psc_trend = A_all @ coef
    psc_phase = _wrap(phase_wrapped - psc_trend)

    psp_az = np.asarray(selection.psp_az[z], dtype=float)
    psp_r = np.asarray(selection.psp_r[z], dtype=float)
    psp_trend = np.column_stack([psp_az, psp_r, np.ones(psp_az.size)]) @ coef
    psp_phase = selection.psp_phase[z]
    if psp_phase is not None and np.asarray(psp_phase).size:
        psp_phase = _wrap(np.asarray(psp_phase, dtype=float) - psp_trend)
    else:
        psp_phase = np.zeros((psp_az.size, n_ifg))

    write_raw(project.path(f"_psc_trend_sel{z + 1}.raw"), psc_trend, "float64")
    write_raw(project.path(f"_psp_trend_sel{z + 1}.raw"), psp_trend, "float64")
    return DetrendResult(
        param=coef, psc_phase=psc_phase, psp_phase=psp_phase,
        psc_trend=psc_trend, psp_trend=psp_trend,
    )


# ---------------------------------------------------------------------------
# group 5/6: densification
# ---------------------------------------------------------------------------


def run_densification(
    project: Project,
    selection: SelectionResult,
    net: NetworkResult,
    *,
    z: int = 0,
    psc_phase: np.ndarray | None = None,
    psp_phase: np.ndarray | None = None,
    psc_h2ph: np.ndarray | None = None,
    psp_h2ph: np.ndarray | None = None,
) -> DensificationResult:
    """Estimate PSP against the nearest PSC and apply the consistency check.

    ``psc_phase``/``psp_phase`` must belong to the same frame as ``net.acheck``
    (i.e. pass the detrended phases when the network was built on detrended data).
    """
    p = project.params
    cand = selection.psc[z]
    if not cand:
        raise ValueError(f"selection {z + 1} contains no PSC")
    net_az = np.array([c.az for c in cand], dtype=float)
    net_r = np.array([c.r for c in cand], dtype=float)
    net_phase = (
        np.stack([c.phase for c in cand]) if psc_phase is None else np.asarray(psc_phase, dtype=float)
    )
    net_h2ph = (
        np.stack([c.h2ph for c in cand]) if psc_h2ph is None else np.asarray(psc_h2ph, dtype=float)
    )
    net_acheck = net.acheck  # NaN for points outside the network

    psp_az = np.asarray(selection.psp_az[z], dtype=float)
    psp_r = np.asarray(selection.psp_r[z], dtype=float)
    psp_phases = selection.psp_phase[z] if psp_phase is None else psp_phase
    psp_h2phs = selection.psp_h2ph[z] if psp_h2ph is None else psp_h2ph
    psp_disp = selection.psp_amp_disp[z]
    if psp_phases is None or np.asarray(psp_phases).size == 0:
        raise ValueError("PSP observations are not available (ps_eval_method != 'psp'?)")

    return densify_points(
        net_az, net_r, net_phase, net_h2ph, net_acheck,
        psp_az, psp_r, np.asarray(psp_phases, dtype=float), np.asarray(psp_h2phs, dtype=float),
        np.asarray(psp_disp, dtype=float),
        Btemp=project.stack.Btemp, Bdop=project.stack.Bdop,
        std_param=p.floats("std_param"), model=int(p.ints("final_model", 1)[0]),
        m2ph=project.m2ph, sat_vel=project.vel,
        nest=p.int("Nest", 3), dens_check=p.str("dens_check", "3of3"),
        max_arc_length=p.float("max_arc_length", 1000.0),
        az_spacing=project.az_spacing, r_spacing=project.r_spacing,
        method=p.str("ps_method", "ils"), sig2_est=net.sig2_est,
        breakpoint=None if p.is_empty("breakpoint") else p.int("breakpoint"),
        breakpoint2=None if p.is_empty("breakpoint2") else p.int("breakpoint2"),
    )


# ---------------------------------------------------------------------------
# group 4: atmosphere / non-linear deformation separation
# ---------------------------------------------------------------------------


def run_filtering(project: Project, net: NetworkResult, *, z: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Split PSC residual phase into non-linear deformation and atmosphere.

    Returns ``(non_linear, atmosphere)`` with the DePSI master-atmosphere term
    (``m2ph * param(2)``) added back to the atmosphere, and writes
    ``_psc_atmo_filt_sel*.raw`` (rows = PSC, columns = interferograms).
    """
    p = project.params
    phase_res = net.psc.phase_res
    master_atmo = None
    if p.ints("final_model", 1)[0] != 1 and net.psc.param.shape[1] > 1:
        master_atmo = project.m2ph * net.psc.param[:, 1]
    non_linear, atmosphere = separate_atmosphere(
        phase_res, project.stack.Btemp,
        length=p.float("ts_atmo_filter_length", 3.0),
        kind=p.str("ts_atmo_filter", "gaussian"),
        master_atmo=master_atmo,
    )
    write_raw(project.path(f"_psc_atmo_filt_sel{z + 1}.raw"), atmosphere, "float64")
    return non_linear, atmosphere


# ---------------------------------------------------------------------------
# group 7/8: consistency, noise filtering, results
# ---------------------------------------------------------------------------


def run_stc(
    project: Project,
    deformation: np.ndarray,
    az: np.ndarray,
    r: np.ndarray,
    *,
    to_mm: bool = False,
) -> np.ndarray:
    """Spatio-temporal consistency of the deformation time series (``ps_spatio_temporal_consistency.m``).

    ``deformation`` is in metres (the ``defo`` block of ``_ps_results_*``), which
    is what the MATLAB implementation reads; the returned value therefore has the
    unit of ``deformation`` (metres by default, millimetres with ``to_mm=True``).
    """
    lo, hi = (float(v) for v in project.params.floats("stc_min_max", [30.0, 100.0]))
    return spatiotemporal_consistency(
        deformation, az, r,
        min_dist=lo, max_dist=hi,
        az_spacing=project.az_spacing, r_spacing=project.r_spacing,
        to_mm=to_mm,
    )


def run_noise_filter(project: Project, defo: np.ndarray) -> np.ndarray:
    """Temporal smoothing of the deformation time series (step 8)."""
    p = project.params
    length = p.float("ts_noise_filter_length", 10.0 / 12.0)
    return filter_deformation(
        defo, project.stack.Btemp,
        length=length, kind=p.str("ts_noise_filter", "gaussian"),
    )


def _ps_results_layout(
    az, r, param, ens_coh, ens_coh_local, sig2hat, phase_unw, defo, atmo,
) -> np.ndarray:
    """Layout read by ``ps_noise_filter``/``ps_output``.

    ``[order az r param(Npar_max) ens_coh ens_coh_local sig2hat phase_unw defo atmo]``
    with one row per PS (``Nifgs`` columns per time-series block).
    """
    n = az.size
    n_ifg = defo.shape[1]
    data = np.empty((n, 6 + NPAR_MAX + 3 * n_ifg))
    data[:, 0] = np.arange(1, n + 1)
    data[:, 1] = az
    data[:, 2] = r
    data[:, 3 : 3 + NPAR_MAX] = param
    data[:, 3 + NPAR_MAX] = ens_coh
    data[:, 4 + NPAR_MAX] = ens_coh_local
    data[:, 5 + NPAR_MAX] = sig2hat
    off = 6 + NPAR_MAX
    data[:, off : off + n_ifg] = phase_unw
    off += n_ifg
    data[:, off : off + n_ifg] = defo
    off += n_ifg
    data[:, off : off + n_ifg] = atmo
    return data


def run_output(
    project: Project,
    selection: SelectionResult,
    net: NetworkResult,
    *,
    z: int = 0,
    densification: DensificationResult | None = None,
) -> dict[str, object]:
    """Write ``_psc_results_sel*``, ``_ps_results_atmo_sel*``, ``_ps_results_defo_sel*``.

    ``_ps_results_*`` is the input of ``ps_noise_filter`` (group 8) and carries
    ``[order az r param(Npar_max) ens_coh ens_coh_local sig2hat phase_unw defo atmo]``.
    """
    cand = selection.psc[z]
    az = np.array([c.az for c in cand], dtype=float)
    r = np.array([c.r for c in cand], dtype=float)
    acheck = net.unwrap.acheck
    phase_unw = net.unwrap.phase_unw
    _, atmosphere = run_filtering(project, net, z=z)

    # _psc_results_sel: [param(Npar_max) ens_coh sig2hat phase_res acheck defo]
    psc_results = np.column_stack(
        [net.psc.param, net.psc.ens_coh, net.psc.sig2hat, net.psc.phase_res, acheck, net.psc.defo]
    )
    write_raw(project.path(f"_psc_results_sel{z + 1}.raw"), psc_results, "float64")

    ens_coh_local = np.abs(np.exp(1j * net.psc.phase_res).mean(axis=1))
    valid = np.isfinite(phase_unw).all(axis=1) & np.isfinite(net.psc.defo).all(axis=1)
    results = _ps_results_layout(
        az[valid], r[valid], net.psc.param[valid], net.psc.ens_coh[valid],
        ens_coh_local[valid], net.psc.sig2hat[valid],
        phase_unw[valid], net.psc.defo[valid], atmosphere[valid],
    )
    write_raw(project.path(f"_ps_results_atmo_sel{z + 1}.raw"), results, "float64")

    geocoded = None
    if project.master_res is not None:
        meta = read_doris_metadata(project.master_res)
        line = az[valid] + project.crop[0] - 1
        pixel = r[valid] + project.crop[2] - 1
        height = net.psc.param[valid, 0] + project.params.float("ref_height", 0.0)
        xyz = lph2xyz(line, pixel, height, meta, orbit_state(meta, line)[1])
        latlonh = xyz_to_latlonh(xyz)
        latlonh[:, :2] = np.rad2deg(latlonh[:, :2])
        geocoded = latlonh
        write_raw(
            project.path(f"_psc_latlon_sel{z + 1}.raw"),
            np.column_stack([az[valid], r[valid], latlonh]),
            "float64",
        )

    stc = run_stc(project, net.psc.defo[valid], az[valid], r[valid])
    stc_filt = run_stc(project, run_noise_filter(project, net.psc.defo[valid]), az[valid], r[valid])

    summary: dict[str, object] = {
        "project_id": project.project_id,
        "selection": z + 1,
        "Nifgs": project.Nifgs,
        "Npsc": int(az.size),
        "Npsp": int(np.asarray(selection.psp_az[z]).size),
        "model": net.model,
        "final_model": net.final_model,
        "Nref": int(net.unwrap.n_ref),
        "groups": int(np.unique(net.unwrap.groups[net.unwrap.groups > 0]).size),
        "adapted_ambiguities": int(net.unwrap.adapted),
        "sig2_est": [float(v) for v in net.sig2_est],
        "sigma0": float(np.nanmedian(net.psc.sig2hat)),
        "stc_median_m": float(np.nanmedian(stc)) if np.isfinite(stc).any() else float("nan"),
        "stc_filt_median_m": float(np.nanmedian(stc_filt)) if np.isfinite(stc_filt).any() else float("nan"),
        "geocoded": bool(geocoded is not None),
    }

    if densification is not None:
        acc = np.asarray(densification.accepted, dtype=bool) & np.isfinite(densification.acheck).all(axis=1)
        idx = np.nonzero(acc)[0]
        if idx.size:
            dens_fit = fit_points(
                densification.phase, densification.h2ph, net.final_model,
                project.stack.Btemp, project.stack.Bdop, project.params.floats("std_param"),
                net.sig2_est, m2ph=project.m2ph, sat_vel=project.vel,
            )
            dens_phase_unw = 2 * np.pi * densification.acheck + np.asarray(densification.phase, dtype=float)
            dens_atmo = np.zeros_like(dens_fit.phase_res)
            dens_layout = _ps_results_layout(
                densification.az[idx], densification.r[idx], dens_fit.param[idx],
                dens_fit.ens_coh[idx], dens_fit.ens_coh[idx], dens_fit.sig2hat[idx],
                dens_phase_unw[idx], dens_fit.defo[idx], dens_atmo[idx],
            )
            # DePSI name for the densified deformation results (read by
            # ps_noise_filter/ps_spatio_temporal_consistency), plus a short alias
            write_raw(project.path(f"_ps_results_defo_sel{z + 1}.raw"), dens_layout, "float64")
            write_raw(project.path(f"_psp_results_sel{z + 1}.raw"), dens_layout, "float64")
            write_point_csv(
                project.path(f"_psp_results_sel{z + 1}.csv"),
                densification.az[idx], densification.r[idx],
                np.column_stack([dens_fit.defo[idx] * 1000.0, densification.votes[idx].astype(float)]),
                columns=[f"defo_{i}" for i in range(project.Nifgs)] + ["votes"],
            )
        summary["densification_candidates"] = int(densification.az.size)
        summary["densification_accepted"] = int(idx.size)
        summary["densification_rate"] = float(idx.size / max(densification.az.size, 1))

    write_resfile(
        project.path("_resfile.txt"),
        {
            "project_id": project.project_id,
            "Npsc": summary["Npsc"],
            "Npsp": summary["Npsp"],
            "model": summary["model"],
            "final_model": summary["final_model"],
            "groups": summary["groups"],
            "adapted_ambiguities": summary["adapted_ambiguities"],
            "stc_median_m": summary["stc_median_m"],
        },
    )
    csv_values = [net.psc.defo[valid] * 1000.0, stc, stc_filt]
    csv_columns = [f"defo_{i}" for i in range(project.Nifgs)] + ["stc_m", "stc_filt_m"]
    if geocoded is not None:
        csv_values.append(geocoded)
        csv_columns += ["lat_deg", "lon_deg", "height_m"]
    write_point_csv(
        project.path(f"_psc_results_sel{z + 1}.csv"),
        az[valid], r[valid],
        np.column_stack(csv_values),
        columns=csv_columns,
    )
    write_point_stack(
        project.path(f"_psc_timeseries_sel{z + 1}.raw"),
        az[valid], r[valid], net.psc.defo[valid],
    )
    return summary


# ---------------------------------------------------------------------------
# orchestrator
# ---------------------------------------------------------------------------


def parse_groups(value) -> tuple[int, ...]:
    """Accept ``"[1:8]"``/``"3,4,6"``/``[3,4]`` and return a tuple of ints."""
    if value is None:
        return tuple(range(1, 9))
    if isinstance(value, (list, tuple, np.ndarray)):
        return tuple(int(v) for v in np.asarray(value).ravel())
    text = str(value).strip().strip("[]")
    if ":" in text or "," in text:
        parts = text.replace(":", ",").split(",")
        nums = [int(float(p)) for p in parts if p.strip()]
        if len(nums) == 3 and ":" in text:
            return tuple(range(nums[0], nums[1] + 1, nums[2] if nums[2] else 1))
        if len(nums) == 2 and ":" in text:
            return tuple(range(nums[0], nums[1] + 1))
        return tuple(nums)
    return (int(float(text)),)


def run_pipeline(
    param_file,
    *,
    workdir=None,
    groups=None,
    progress=None,
) -> dict[str, object]:
    """Run the requested DePSI groups and return a summary of the run.

    ``groups`` follows ``processing_groups`` (1..8); groups 5 and 6 both mean
    "densification" (the MATLAB driver runs it in group 5 and again in group 6
    for the deformation model, which is disabled by default).
    """
    groups = parse_groups(groups)
    started = time.perf_counter()
    project = load_project(param_file, workdir=workdir)
    summary: dict[str, object] = {
        "project_id": project.project_id,
        "groups": list(groups),
        "Nlines": project.Nlines,
        "Npixels": project.Npixels,
        "Nifgs": project.Nifgs,
        "Nslc": project.Nslc,
    }

    def note(stage: str) -> None:
        if progress is not None:
            progress(stage, time.perf_counter() - started)

    selection = None
    if 2 in groups or any(g in groups for g in (3, 4, 5, 6, 7, 8)):
        selection = run_selection(project)
        summary["Npsc"] = [int(v) for v in np.atleast_1d(selection.Npsc)]
        summary["Npsp"] = [int(v) for v in np.atleast_1d(selection.Npsp)]
        note("selection")

    if selection is None:
        if any(g in groups for g in (2, 3, 4, 5, 6, 7, 8)):
            raise ValueError("a processing run needs group 2 (selection) at minimum")
        summary["wall_time_s"] = round(time.perf_counter() - started, 3)
        write_summary(project.path("_summary.json"), summary)
        return summary

    net = None
    detrend = None
    if 3 in groups:
        net = run_network(project, selection)
        summary["Npsc_network"] = int(np.isfinite(net.acheck).all(axis=1).sum())
        summary["sig2_est"] = [float(v) for v in net.sig2_est]
        summary["adapted_ambiguities"] = int(net.unwrap.adapted)
        note("network")
        if project.params.str("detrend_method", "yes").lower() == "yes":
            detrend = run_detrend(project, selection, net)
            net_detrended = run_network(project, selection, phase=detrend.psc_phase)
            summary["adapted_ambiguities_after_detrend"] = int(net_detrended.unwrap.adapted)
            # The per-IFG plane need not be compatible with the height-to-phase
            # scaling, so on topography-dominated data the second network can lose
            # integer consistency.  Only accept the detrended frame when the
            # network stayed consistent and the fit did not degrade.
            cap = int(net_detrended.unwrap.arcs.shape[0] * project.Nifgs)
            sigma_before = float(np.nanmedian(net.psc.sig2hat))
            sigma_after = float(np.nanmedian(net_detrended.psc.sig2hat))
            accepted = (
                net_detrended.unwrap.adapted < 0.5 * max(cap, 1)
                and np.isfinite(sigma_after)
                and sigma_after <= 5.0 * sigma_before
            )
            summary["detrend_accepted"] = accepted
            summary["sigma0_after_detrend"] = sigma_after
            if accepted:
                net = net_detrended
            else:
                detrend = None
            note("detrend")

    if 4 in groups:
        if net is None:
            net = run_network(project, selection)
        non_linear, atmosphere = run_filtering(project, net)
        summary["atmo_std_rad"] = float(np.nanstd(atmosphere))
        summary["non_linear_std_rad"] = float(np.nanstd(non_linear))
        note("filtering")

    densification = None
    if 5 in groups or 6 in groups:
        if net is None:
            net = run_network(project, selection)
        densification = run_densification(
            project, selection, net,
            psc_phase=None if detrend is None else detrend.psc_phase,
            psp_phase=None if detrend is None else detrend.psp_phase,
        )
        summary["densification_candidates"] = int(densification.az.size)
        summary["densification_accepted"] = int(densification.accepted.sum())
        note("densification")

    if 7 in groups or 8 in groups:
        if net is None:
            net = run_network(project, selection)
        out = run_output(project, selection, net, densification=densification)
        summary.update({k: v for k, v in out.items() if k not in summary})
        note("output")

    summary["wall_time_s"] = round(time.perf_counter() - started, 3)
    write_summary(project.path("_summary.json"), summary)
    return summary
