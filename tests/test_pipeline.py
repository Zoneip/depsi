"""End-to-end pipeline tests (groups 2..8) on a small synthetic project."""

from __future__ import annotations

import json

import numpy as np
import pytest

from depsipy.cli import main
from depsipy.io.raw import RawFile
from depsipy.pipeline import (
    load_project,
    run_densification,
    run_detrend,
    run_network,
    run_noise_filter,
    run_output,
    run_pipeline,
    run_selection,
    run_stc,
)
from depsipy.simulate import generate_project


@pytest.fixture(scope="module")
def project_dir(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("pipeline")
    sim = generate_project(
        workdir,
        Nlines=60,
        Npixels=60,
        Nslc=15,
        seed=2,
        ps_fraction=0.15,
        atmo_amp=0.2,
        defo_rate=0.03,
        gridsize=5.0,
    )
    return sim


@pytest.fixture(scope="module")
def truth(project_dir):
    return dict(np.load(project_dir.truth_file))


def _truth_phase(project, selection, truth, z=0, detrend=None):
    """Unwrapped phase of every PSC in the interferometric frame."""
    cand = selection.psc[z]
    az = np.array([c.az for c in cand])
    r = np.array([c.r for c in cand])
    slaves = truth["slave_index"]
    m2ph = float(truth["m2ph"][0])
    dem = truth["dem"]
    h2ph_scene = truth["h2ph_scene"][slaves]
    phase = (m2ph * h2ph_scene[:, None] * dem[az - 1, r - 1][None, :] +
             m2ph * truth["defo"][slaves][:, az - 1, r - 1]).T
    if detrend is not None:
        phase = phase - detrend.psc_trend
    return phase


def test_step3_recovers_the_true_phase(project_dir, truth):
    project = load_project(project_dir.param_file)
    selection = run_selection(project)
    net = run_network(project, selection)
    ref = int(net.unwrap.ref_nodes[0])
    true_phase = _truth_phase(project, selection, truth)
    est = net.unwrap.phase_unw
    ok = np.isfinite(est).all(axis=1)
    resid = est[ok] - (true_phase - true_phase[ref])[ok]
    assert ok.sum() >= 0.9 * est.shape[0]
    assert np.median(np.abs(resid)) < 1.0
    assert net.unwrap.adapted == 0


def test_detrend_keeps_the_phase_consistent(project_dir, truth):
    project = load_project(project_dir.param_file)
    selection = run_selection(project)
    net = run_network(project, selection)
    detrend = run_detrend(project, selection, net)
    net2 = run_network(project, selection, phase=detrend.psc_phase)

    true_phase = _truth_phase(project, selection, truth, detrend=detrend)
    ref = int(net2.unwrap.ref_nodes[0])
    est = net2.unwrap.phase_unw
    ok = np.isfinite(est).all(axis=1)
    resid = est[ok] - (true_phase - true_phase[ref])[ok]
    assert np.median(np.abs(resid)) < 1.0
    assert net2.unwrap.adapted == 0
    # the plane fit must not collapse onto a degenerate subset
    assert np.nanmax(np.abs(detrend.param[:2])) < 20.0


def test_densification_enriches_true_ps(project_dir, truth):
    project = load_project(project_dir.param_file)
    selection = run_selection(project)
    net = run_network(project, selection)
    detrend = run_detrend(project, selection, net)
    net = run_network(project, selection, phase=detrend.psc_phase)
    dens = run_densification(
        project, selection, net,
        psc_phase=detrend.psc_phase, psp_phase=detrend.psp_phase,
    )
    accepted = np.asarray(dens.accepted, dtype=bool)
    assert accepted.sum() > 0
    az, r = dens.az[accepted], dens.r[accepted]
    ps_mask = truth["ps_mask"]
    # the candidates are mostly distributed scatterers; accepted ones must be
    # strongly biased towards the true persistent scatterers
    assert ps_mask[az - 1, r - 1].mean() > ps_mask.mean() + 0.2
    assert np.median(dens.varfac[accepted]) < 3.0


def test_noise_filter_and_stc(project_dir):
    project = load_project(project_dir.param_file)
    selection = run_selection(project)
    net = run_network(project, selection)
    filtered = run_noise_filter(project, net.psc.defo)
    assert filtered.shape == net.psc.defo.shape
    ok = np.isfinite(net.psc.defo).all(axis=1)
    assert np.isfinite(filtered[ok]).all()
    # smoothing must reduce the temporal roughness of the series
    rough_before = np.abs(np.diff(net.psc.defo[ok], axis=1)).mean()
    rough_after = np.abs(np.diff(filtered[ok], axis=1)).mean()
    assert rough_after < rough_before

    project.params.values["stc_min_max"] = np.array([1.0, 60.0])
    az = np.array([c.az for c in selection.psc[0]], float)
    r = np.array([c.r for c in selection.psc[0]], float)
    stc = run_stc(project, net.psc.defo, az, r)
    assert np.isfinite(stc).mean() > 0.5


def test_run_pipeline_writes_the_depsi_products(project_dir):
    workdir = project_dir.workdir
    summary = run_pipeline(project_dir.param_file, groups=(2, 3, 4, 6, 7, 8))
    for suffix in (
        "_psc_2orig_sel1.raw",
        "_psp_2orig_sel1.raw",
        "_psc_results_sel1.raw",
        "_ps_results_atmo_sel1.raw",
        "_psp_results_sel1.raw",
        "_psc_atmo_filt_sel1.raw",
        "_psc_timeseries_sel1.raw",
        "_summary.json",
        "_resfile.txt",
    ):
        assert (workdir / f"sim{suffix}").exists(), suffix
    assert summary["Npsc"] == [int(summary["Npsc"][0])]
    assert summary["densification_accepted"] > 0

    # _psc_results_sel1 carries [param(Npar_max) ens_coh sig2hat phase_res acheck defo]
    n_ifg = summary["Nifgs"]
    path = workdir / "sim_psc_results_sel1.raw"
    data = np.fromfile(path, dtype=np.float64)
    assert data.size == summary["Npsc"][0] * (2 + 10 + 3 * n_ifg)

    # _psp_results_sel1 is the layout ps_noise_filter reads
    path = workdir / "sim_psp_results_sel1.raw"
    data = np.fromfile(path, dtype=np.float64)
    assert data.size == summary["densification_accepted"] * (6 + 10 + 3 * n_ifg)

    with open(workdir / "sim_summary.json", encoding="utf-8") as fh:
        assert json.load(fh)["Npsc"] == summary["Npsc"]


def _write_geocoding_res(path, *, lines: int, samples: int) -> Path:
    """A Doris ``.res`` whose orbit covers the synthetic project's line numbers."""
    a = 6378137.0
    radius = a + 700e3
    t = np.linspace(0.0, 1.0, 41)
    orbit = np.column_stack([t, np.full_like(t, radius), 7000.0 * t, np.zeros_like(t)])
    head = [
        " MASTER RESULTFILE:  master.res",
        "Radar_wavelength (m):                           0.0555",
        "First_pixel_azimuth_time (UTC):                 01-JAN-2020 00:00:00.000",
        "Pulse_Repetition_Frequency (computed, Hz):      1000.0",
        "Range_time_to_first_pixel (2way) (ms):          5.4",
        "Range_sampling_rate (computed, MHz):            30.0",
        "Scene_centre_latitude:                          2.86",
        "Scene_centre_longitude:                         0.0",
        "Xtrack_f_DC_constant (Hz, early edge):          0.0",
        "Xtrack_f_DC_linear (Hz/s, early edge):          0.0",
        "Xtrack_f_DC_quadratic (Hz/s/s, early edge):     0.0",
        f"Number_of_lines_original:                       {lines}",
        f"Number_of_pixels_original:                      {samples}",
        "First_line (w.r.t. original_image):             1",
        f"Last_line (w.r.t. original_image):              {lines}",
        "First_pixel (w.r.t. original_image):            1",
        f"Last_pixel (w.r.t. original_image):             {samples}",
        f"NUMBER_OF_DATAPOINTS:                           {orbit.shape[0]}",
    ]
    rows = [f"{row[0]:.6f} {row[1]:.4f} {row[2]:.4f} {row[3]:.4f}" for row in orbit]
    path.write_text("\n".join(head + rows) + "\n", encoding="utf-8")
    return path


def test_output_geocodes_when_master_res_is_given(project_dir, truth):
    project = load_project(project_dir.param_file)
    selection = run_selection(project)
    net = run_network(project, selection)
    res = _write_geocoding_res(
        project_dir.workdir / "master.res", lines=project.Nlines, samples=project.Npixels
    )
    project.master_res = str(res)
    summary = run_output(project, selection, net)
    assert summary["geocoded"] is True

    latlon_path = project_dir.workdir / "sim_psc_latlon_sel1.raw"
    data = np.fromfile(latlon_path, dtype=np.float64).reshape(-1, 5)
    assert data.shape[0] > 0
    assert np.isfinite(data).all()
    lat, lon, height = data[:, 2], data[:, 3], data[:, 4]
    assert np.all((lat > -90) & (lat < 90))
    assert np.all((lon >= -180) & (lon <= 180))
    # the recovered height must agree with the topographic parameter used as input
    param = net.psc.param[:, 0][np.isfinite(net.unwrap.phase_unw).all(axis=1)]
    assert np.max(np.abs(height - param)) < 0.01

    csv = np.genfromtxt(
        project_dir.workdir / "sim_psc_results_sel1.csv", delimiter=",", names=True
    )
    assert {"lat_deg", "lon_deg", "height_m"} <= set(csv.dtype.names)


def test_cli_runs_and_reports_json(project_dir, capsys):
    workdir = project_dir.workdir
    rc = main([
        "run", str(project_dir.param_file),
        "--workdir", str(workdir),
        "--groups", "2,3",
        "--json", "--quiet",
    ])
    out = capsys.readouterr().out
    assert rc == 0
    payload = json.loads(out)
    assert payload["project_id"] == "sim"
    assert payload["groups"] == [2, 3]
    assert payload["Npsc"][0] > 0
