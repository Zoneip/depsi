"""Synthetic DePSI project generator.

Creates a complete, self-contained data set in the format the MATLAB toolbox
expects (``*_slc.raw``, ``*_ifg.raw``, ``*_h2ph.raw`` + ``input_file.txt`` +
``param_file.txt``) together with the ground truth, so that the whole processing
chain can be validated against known parameters.

This replaces the 1050-line ``depsi_sim.m`` (which writes a DePSI-specific
mat-file workspace) with a compact generator whose output can be fed straight to
``depsipy run`` or to the original MATLAB code.

Physics
-------
For every scene ``i`` with temporal baseline ``Btemp_i`` and height-to-phase
factor ``h2ph_i`` the interferometric phase (master removed) is

    phase_i(p) = m2ph * h2ph_i * h(p)            topography
               + m2ph * d_i(p)                   deformation
               + a_i(p)                          atmospheric screen
               + n_i(p)                          noise

with ``m2ph = -4*pi/lambda`` and ``d`` in metres.  ``h2ph`` is written as a
separate raster because that is what DePSI reads; keeping it consistent with the
generated phase is what makes the test meaningful.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["SyntheticScene", "SyntheticProject", "generate_project"]


@dataclass
class SyntheticScene:
    date: str
    orbit: str
    Bperp: float
    Btemp: float
    Bdop: float


@dataclass
class SyntheticProject:
    """Paths and truth of a generated project."""

    workdir: Path
    param_file: Path
    input_file: Path
    truth_file: Path
    scenes: list[SyntheticScene]
    master_index: int
    Nlines: int
    Npixels: int
    h2ph: np.ndarray
    truth: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# field generators
# ---------------------------------------------------------------------------


def _fractal_field(rng, shape, beta=2.5, scale=1.0):
    """Spatially correlated field with a 1/f^beta power spectrum (fft filtering)."""
    n0, n1 = shape
    ky = np.fft.fftfreq(n0)[:, None]
    kx = np.fft.fftfreq(n1)[None, :]
    k = np.sqrt(ky**2 + kx**2)
    k[0, 0] = 1.0
    amp = k ** (-beta / 2.0)
    amp[0, 0] = 0.0
    spec = (rng.normal(size=shape) + 1j * rng.normal(size=shape)) * amp
    f = np.real(np.fft.ifft2(spec))
    f -= f.mean()
    s = f.std()
    return f / s * scale if s > 0 else f


def _smooth_dem(rng, shape, scale=60.0):
    """A smooth synthetic DEM in metres."""
    n0, n1 = shape
    y = np.linspace(0, 1, n0)[:, None]
    x = np.linspace(0, 1, n1)[None, :]
    dem = (
        20.0 * np.sin(2 * np.pi * x) * np.cos(3 * np.pi * y)
        + 10.0 * np.sin(5 * np.pi * x + 1.0)
        + _fractal_field(rng, shape, beta=3.0, scale=6.0)
    )
    return dem * (scale / 20.0)


def _deformation(rng, shape, Btemp, rate, seasonal, event):
    """Spatiotemporal deformation, in metres, shape (Nifgs, Nlines, Npixels)."""
    n = Btemp.size
    n0, n1 = shape
    y = np.linspace(0, 1, n0)[:, None]
    x = np.linspace(0, 1, n1)[None, :]
    # a subsidence bowl centred off-centre plus a smooth regional gradient
    bowl = np.exp(-(((x - 0.4) ** 2 + (y - 0.55) ** 2) / 0.05))
    spatial = (0.6 * bowl - 0.4 * x + 0.2 * y) * rate  # m
    out = np.empty((n, n0, n1))
    phase_year = np.linspace(-2.0, 2.0, n)  # years relative to the master epoch
    for i, t in enumerate(Btemp):
        out[i] = spatial * t if t != 0 else 0.0
    if seasonal:
        # annual signal, 1 cm amplitude, phase locked to the calendar
        out += seasonal * np.cos(2 * np.pi * (phase_year[:, None, None]))[..., None, :] * 0
        out += seasonal * np.sin(2 * np.pi * np.linspace(0, 4, n))[:, None, None] * spatial / rate
    if event:
        step_index, step_size = event
        out[step_index:] += step_size * bowl
    out = out - out[[int(np.argmin(np.abs(Btemp)))]]  # master epoch is the zero
    return out


# ---------------------------------------------------------------------------
# project generation
# ---------------------------------------------------------------------------


def generate_project(
    workdir,
    *,
    Nlines: int = 120,
    Npixels: int = 120,
    Nslc: int = 12,
    seed: int = 0,
    az_spacing: float = 1.8,
    r_spacing: float = 0.62,
    sensor: str = "spacety",
    Bperp: float | None = None,
    Btemp_days: float = 24.0,
    m2ph: float = -4 * np.pi / 0.055517120823757024,
    defo_rate: float = 0.05,
    seasonal: float = 0.0,
    event: tuple[int, float] | None = None,
    atmo_amp: float = 0.6,
    noise_amp: float = 0.25,
    ps_fraction: float = 0.08,
    amplitude: float = 1000.0,
    nlooks_phase: float = 4.0,
    write_param_file: bool = True,
    processing_groups: str = "[1:8]",
    gridsize: float | None = None,
) -> SyntheticProject:
    """Generate a synthetic DePSI project on disk.

    Parameters
    ----------
    workdir
        Directory to create (must be empty or non-existent).
    Nlines, Npixels
        Raster size.
    Nslc
        Number of SLCs (one of them is the master).
    ps_fraction
        Fraction of pixels that are "persistent": they keep a stable amplitude
        in every scene.  The rest are distributed scatterers whose amplitude
        varies strongly from scene to scene (what the PSC selection has to
        reject).
    noise_amp
        Phase noise in radians (standard deviation).
    atmo_amp
        Standard deviation of the atmospheric screen in radians.
    """
    rng = np.random.default_rng(seed)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    shape = (Nlines, Npixels)
    Nifgs = Nslc - 1

    # -- acquisition geometry ------------------------------------------------
    if Bperp is None:
        Bperp = np.round(np.linspace(-600.0, 600.0, Nslc) + rng.normal(0, 20, Nslc), 1)
    else:
        Bperp = np.asarray(Bperp, dtype=float)

    # temporal baselines in years, master at index Nslc//2
    master = Nslc // 2
    days = (np.arange(Nslc) - master) * Btemp_days
    Btemp_days_arr = days.astype(float)
    Btemp = Btemp_days_arr / 365.25
    Bdop = np.round(rng.normal(0, 10, Nslc), 1)

    # height-to-phase factors (1/m), one per scene
    R = 800e3
    inc = np.deg2rad(23.0)
    h2ph_scene = -2.0 * Bperp / (0.055517120823757024 * R * np.sin(inc))

    # -- ground truth --------------------------------------------------------
    dem = _smooth_dem(rng, shape)
    defo = _deformation(rng, shape, Btemp, defo_rate, seasonal, event)

    ps_mask = np.zeros(shape, dtype=bool)
    n_ps = int(ps_fraction * Nlines * Npixels)
    idx = rng.choice(Nlines * Npixels, size=n_ps, replace=False)
    ps_mask.flat[idx] = True
    # a PS "cluster" so that the network has a dense core like real data
    ps_mask[Nlines // 6 : Nlines // 3, Npixels // 3 : Npixels // 2] |= True

    amp_ps = amplitude * (1.0 + 0.1 * rng.normal(size=shape))
    amp_ds = amplitude * np.abs(rng.normal(size=shape) + 0.6)
    amp_field = np.where(ps_mask, amp_ps, amp_ds)

    # per-pixel, per-scene amplitude modulation of distributed scatterers
    ds_modulation = np.abs(rng.normal(1.0, 0.45, size=(Nslc,) + shape))
    atmo_scene = np.array([_fractal_field(rng, shape, beta=2.7, scale=atmo_amp) for _ in range(Nslc)])

    # -- write rasters -------------------------------------------------------
    filenames_slc, filenames_ifg, filenames_h2ph = [], [], []
    dates, orbits = [], []
    from .io.raw import write_raw

    slc_phases = np.empty((Nslc,) + shape)
    for i in range(Nslc):
        phase = m2ph * h2ph_scene[i] * dem + m2ph * defo[i] + atmo_scene[i]
        slc_phases[i] = phase
        amp = np.where(ps_mask, amp_field, amp_field * ds_modulation[i])
        noise = rng.normal(0, noise_amp / np.sqrt(nlooks_phase) * (1.0 + 3.0 * ~ps_mask), shape)
        slc = (amp * np.exp(1j * (phase + noise))).astype(np.complex64)
        p = workdir / f"slc_{i:03d}.raw"
        write_raw(p, slc, "cpxfloat32")
        filenames_slc.append(p.name)
        dates.append(f"{(i * int(Btemp_days)):05d}")
        orbits.append(str(10000 + i))

    for i in range(Nifgs):
        phase = m2ph * h2ph_scene[i] * dem + m2ph * defo[i] + atmo_scene[i]
        amp = amp_field * np.where(ps_mask, 1.0, ds_modulation[i] * ds_modulation[master])
        # interferometric noise grows with the number of looks lost
        noise = rng.normal(
            0, noise_amp * np.sqrt(2) / np.sqrt(nlooks_phase) * (1.0 + 3.0 * ~ps_mask), shape
        )
        ifg = (amp * np.exp(1j * (phase + noise))).astype(np.complex64)
        p = workdir / f"ifg_{i:03d}.raw"
        write_raw(p, ifg, "cpxfloat32")
        filenames_ifg.append(p.name)

        hp = (h2ph_scene[i] * np.ones(shape)).astype(np.float32)
        ph = workdir / f"h2ph_{i:03d}.raw"
        write_raw(ph, hp, "float32")
        filenames_h2ph.append(ph.name)

    # -- input file ----------------------------------------------------------
    crop = (1, Nlines, 1, Npixels)
    lines = []
    for i in range(Nifgs):
        # DePSI sorts on Btemp and expects the master as the last line
        lines.append(
            f"{orbits[i]} {dates[i]} {filenames_slc[i]} {filenames_ifg[i]} "
            f"{filenames_h2ph[i]} {Path(filenames_ifg[i]).stem} "
            f"{Bperp[i]:.1f} {Btemp[i]:.9f} {Bdop[i]:.1f} "
            f"{crop[0]} {crop[1]} {crop[2]} {crop[3]}"
        )
    lines.append(
        f"{orbits[master]} {dates[master]} {filenames_slc[master]} dummy dummy dummy "
        f"0 0 0 {crop[0]} {crop[1]} {crop[2]} {crop[3]}"
    )
    input_file = workdir / "input.txt"
    input_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # -- truth ---------------------------------------------------------------
    truth = dict(
        dem=dem,
        defo=defo,  # (Nifgs, Nlines, Npixels) metres
        ps_mask=ps_mask,
        h2ph_scene=h2ph_scene,
        Btemp=Btemp,
        Bperp=np.asarray(Bperp, dtype=float),
        Btemp_days=Btemp_days_arr,
        master_index=master,
        m2ph=m2ph,
    )
    truth_file = workdir / "truth.npz"
    np.savez_compressed(
        truth_file,
        dem=dem,
        defo=defo,
        ps_mask=ps_mask,
        h2ph_scene=h2ph_scene,
        Btemp=Btemp,
        Bperp=np.asarray(Bperp, dtype=float),
        Btemp_days=Btemp_days_arr,
        m2ph=np.array([m2ph]),
        master_index=np.array([master]),
    )

    param_file = workdir / "param_file.txt"
    if write_param_file:
        param_file.write_text(
            _param_text(
                workdir=workdir,
                input_file=input_file,
                sensor=sensor,
                az_spacing=az_spacing,
                r_spacing=r_spacing,
                processing_groups=processing_groups,
                master_res=None,
                gridsize=gridsize if gridsize is not None else 10 * az_spacing,
            ),
            encoding="utf-8",
        )

    scenes = [
        SyntheticScene(
            date=dates[i],
            orbit=orbits[i],
            Bperp=float(Bperp[i]),
            Btemp=float(Btemp[i]),
            Bdop=float(Bdop[i]),
        )
        for i in range(Nslc)
    ]

    return SyntheticProject(
        workdir=workdir,
        param_file=param_file,
        input_file=input_file,
        truth_file=truth_file,
        scenes=scenes,
        master_index=master,
        Nlines=Nlines,
        Npixels=Npixels,
        h2ph=h2ph_scene,
        truth=truth,
    )


def _param_text(*, workdir, input_file, sensor, az_spacing, r_spacing, processing_groups, master_res, gridsize):
    return f"""Input_file template for depsipy
% General parameters
max_mem_buffer = 50e7
visible_plots = 'n'
detail_plots = 'n'
processing_groups = {processing_groups}
run_mode = 'normal'

% Project parameters
project_id = 'sim'
input_file = '{input_file}'
sensor = '{sensor}'
orbit = 'asc'
crop = []
az_spacing = {az_spacing}
r_spacing = {r_spacing}
slc_selection_input = []
ifg_selection_input = []
ref_cn = []
Ncv = 25
ps_method = 'boot'
psc_model = [1]
ps_model = [1]
final_model = [2]
breakpoint = []
breakpoint2 = []
ens_coh_threshold = 0.9
varfac_threshold = 3
detrend_method = 'yes'
output_format = 1
stc_min_max = [30,100]
do_apriori_sidelobe_mask = 'yes'
do_aposteriori_sidelobe_mask = 'no'
master_res = {("'" + str(master_res) + "'") if master_res else "[]"}
ref_height = 0

% PSC parameters
amplitude_calibration = 'no'
psc_selection_method = 'threshold'
psc_selection_gridsize = {gridsize}
psc_threshold = 0.3
max_arc_length = 1000
network_method = 'delaunay'
Ncon = 8
Nparts = 8
Npsc_selections = 1
filename_water_mask = []
gamma_threshold = 0.75
psc_distribution = 'nonuniform'
weighted_unwrap = 'yes'
livetime_threshold = 0.2
peak_tolerance = 0.99

% Ps parameters
psp_selection_method = 'ampdisp'
psp_threshold1 = 0.38
psp_threshold2 = []
ps_eval_method = 'psp'
Namp_disp_bins = 100
Ndens_iterations = 5
densification_flag = 'yes'
ps_area_of_interest = []
dens_method = 'orig'
dens_check = '3of3'
Nest = 3

% Stochastic model
std_param = [30,0.005,1,0.02,0.01,0.01,0.005]
defo_range = 5000
weighting = 'vce'
ts_atmo_filter = 'gaussian'
ts_atmo_filter_length = 3
ts_noise_filter = 'gaussian'
ts_noise_filter_length = 10/12

% Bowl parameters
defo_method = []
xc0 = []
yc0 = []
zc0 = []
r0 = []
r10 = []
epoch = []
"""
