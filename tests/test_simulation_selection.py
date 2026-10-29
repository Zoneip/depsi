"""Tests for the synthetic project generator and the PSC/PSP selection.

These are the end-to-end invariants of milestone M1: the generator must produce a
project that the readers can consume, and the selection must find the injected
persistent scatterers without dragging in distributed ones.
"""

from __future__ import annotations

import numpy as np
import pytest

from depsipy.io.params import load_params
from depsipy.io.project import read_input_file
from depsipy.io.raw import infer_shape, read_pixels, read_raw
from depsipy.simulate import generate_project
from depsipy.steps.selection import grid_arrays, select_psc_psp


@pytest.fixture(scope="module")
def project(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("sim")
    return generate_project(
        workdir,
        Nlines=80,
        Npixels=80,
        Nslc=8,
        seed=7,
        ps_fraction=0.10,
        noise_amp=0.25,
    )


def test_generated_files_are_readable(project):
    st = read_input_file(project.param_file.parent / "input.txt")
    assert st.Nifgs == 7
    assert st.Nslc == 8
    assert len(st.slc) == 8  # sorted scenes + master
    assert st.crop_in.shape == (8, 4)

    # complex SLC, real h2ph, and the width follows from the file size
    lines, width = infer_shape(project.workdir / st.slc[0], 80, "cpxfloat32")
    assert (lines, width) == (80, 80)
    slc = read_raw(project.workdir / st.slc[0], 80, "cpxfloat32")
    assert np.iscomplexobj(slc)
    h2ph = read_raw(project.workdir / st.h2ph[0], 80, "float32")
    assert np.allclose(h2ph, h2ph.flat[0])  # constant per interferogram


def test_truth_is_self_consistent(project):
    truth = np.load(project.truth_file)
    assert truth["dem"].shape == (80, 80)
    # deformation is stored per acquisition (Nslc), the interferograms are Nslc-1
    assert truth["defo"].shape == (8, 80, 80)
    assert 0.02 < truth["ps_mask"].mean() < 0.4
    # deformation must vanish at the master epoch
    Btemp = truth["Btemp"]
    master = int(np.argmin(np.abs(Btemp)))
    assert np.abs(truth["defo"][master]).max() < 1e-6
    # a metre-scale DEM and millimetre-to-centimetre deformation
    assert 1.0 < np.abs(truth["dem"]).max() < 500.0
    assert 1e-4 < np.abs(truth["defo"]).max() < 0.5


def test_param_file_roundtrip(project):
    prm = load_params(project.param_file)
    assert prm.str("sensor") == "spacety"
    assert prm.int("Npsc_selections") == 1
    np.testing.assert_array_equal(prm.ints("psc_model"), [1])
    assert prm.float("psc_selection_gridsize") == pytest.approx(18.0)


def test_grid_arrays_cover_the_raster():
    az, r = grid_arrays(125, 120, 1.8, 0.62, 18.0)
    assert az[0, 0] == 1 and az[-1, 1] == 125
    assert r[0, 0] == 1 and r[-1, 1] == 120
    # borders are contiguous and non-overlapping
    for edges in (az, r):
        np.testing.assert_array_equal(edges[1:, 0], edges[:-1, 1] + 1)


def test_selection_finds_the_injected_ps(project):
    st = read_input_file(project.param_file.parent / "input.txt")
    prm = load_params(project.param_file)
    Nlines = st.crop_final[1] - st.crop_final[0] + 1
    Npixels = st.crop_final[3] - st.crop_final[2] + 1

    res = select_psc_psp(
        slc_files=st.slc_paths(),
        ifg_files=[st.resolve(f) for f in st.ifg],
        h2ph_files=[st.resolve(f) for f in st.h2ph],
        Nlines=Nlines,
        Npixels=Npixels,
        az_spacing=prm.float("az_spacing"),
        r_spacing=prm.float("r_spacing"),
        crop=st.crop_final,
        crop_in=st.crop_in,
        slc_selection=np.arange(1, st.Nslc + 1),
        calfactors=np.ones(st.Nslc),
        gridsize=prm.float("psc_selection_gridsize"),
        psc_threshold=prm.float("psc_threshold"),
        psp_threshold1=prm.float("psp_threshold1"),
        psp_threshold2=None,
        psc_distribution=prm.str("psc_distribution"),
        do_apriori_sidelobe_mask=prm.str("do_apriori_sidelobe_mask"),
        peak_tolerance=0.7,
    )

    truth = np.load(project.truth_file)
    mask = truth["ps_mask"]
    assert res.Npsc[0] > 0

    az = np.array([c.az for c in res.psc[0]])
    r = np.array([c.r for c in res.psc[0]])
    on_ps = mask[az - 1, r - 1]

    # precision: nearly every selected PSC must be an injected persistent scatterer
    assert on_ps.mean() >= 0.95, f"precision {on_ps.mean():.2f}"

    # coverage: every grid cell that contains a "reachable" PS must yield a PSC
    reach = (res.amp_disp < prm.float("psc_threshold")) & mask & ~res.side_lobe
    cells_with = set()
    for a0, a1 in res.grid_az:
        for r0, r1 in res.grid_r:
            if reach[a0 - 1 : a1, r0 - 1 : r1].any():
                cells_with.add((a0, r0))
    cells_got = set()
    for a0, a1 in res.grid_az:
        for r0, r1 in res.grid_r:
            if np.any((az >= a0) & (az <= a1) & (r >= r0) & (r <= r1)):
                cells_got.add((a0, r0))
    coverage = len(cells_got & cells_with) / max(len(cells_with), 1)
    assert coverage >= 0.95, f"coverage {coverage:.2f}"

    # the amplitude dispersion must separate the two populations
    med_ps = np.nanmedian(res.amp_disp[mask])
    med_ds = np.nanmedian(res.amp_disp[~mask])
    assert med_ps < 0.25 < med_ds

    # the side-lobe mask must not eat the persistent scatterers
    assert res.side_lobe[mask].mean() < 0.05


def test_selection_amplitude_maps_are_plausible(project):
    st = read_input_file(project.param_file.parent / "input.txt")
    prm = load_params(project.param_file)
    res = select_psc_psp(
        slc_files=st.slc_paths(),
        ifg_files=[st.resolve(f) for f in st.ifg],
        h2ph_files=[st.resolve(f) for f in st.h2ph],
        Nlines=80,
        Npixels=80,
        az_spacing=prm.float("az_spacing"),
        r_spacing=prm.float("r_spacing"),
        crop=st.crop_final,
        crop_in=st.crop_in,
        slc_selection=np.arange(1, st.Nslc + 1),
        calfactors=np.ones(st.Nslc),
        gridsize=prm.float("psc_selection_gridsize"),
        psc_threshold=prm.float("psc_threshold"),
        psp_threshold1=prm.float("psp_threshold1"),
        psp_threshold2=None,
        psc_distribution="uniform",
        do_apriori_sidelobe_mask="no",
        peak_tolerance=0.7,
    )
    assert res.mean_amp.shape == (80, 80)
    assert np.all(np.isfinite(res.mean_amp))
    assert res.amp_disp.shape == (80, 80)
    # the multi-reflectivity map is a log10 amplitude: ~1e2..1e3 -> 2..3
    assert 1.5 < np.nanmedian(np.log10(res.mean_amp)) < 3.5
