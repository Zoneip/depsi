"""Same-data parity against the MATLAB sources executed by Octave.

The MATLAB toolbox cannot run on this machine (no MATLAB; the ARM node is
x86-only for MATLAB), so ``ps_selection.m`` is executed *unmodified* through
Octave 9.4 + the shims in ``benchmarks/octave/shim``.  The test runs both
implementations on one synthetic project and requires identical PSC/PSP lists.

Skipped when Octave or the local ``reference``/``workspace`` toolbox is absent
(i.e. on a clean checkout or on the cloud node).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from octave_parity import (  # noqa: E402
    compare_points,
    read_points,
    run_octave_selection,
    write_selection_input,
)

from depsipy.pipeline import load_project, run_selection  # noqa: E402
from depsipy.simulate import generate_project  # noqa: E402

TOOLBOX = REPO / "workspace/ps_toolbox_v1.7.8.8_hhf"


@pytest.mark.skipif(shutil.which("octave") is None, reason="Octave is not installed")
@pytest.mark.skipif(
    not (TOOLBOX / "main/ps_selection.m").exists(), reason="MATLAB toolbox is local-only"
)
def test_selection_matches_matlab_on_the_same_project(tmp_path):
    sim = generate_project(
        tmp_path / "sim",
        Nlines=60, Npixels=60, Nslc=15, seed=2,
        ps_fraction=0.15, atmo_amp=0.2, defo_rate=0.03, gridsize=5.0,
    )
    scenario = tmp_path / "octave_input"
    scenario.mkdir()
    write_selection_input(sim.param_file, scenario)
    workdir = tmp_path / "octave_work"
    workdir.mkdir()

    result = run_octave_selection(TOOLBOX, scenario, workdir, octave=shutil.which("octave"))
    assert result["returncode"] == 0, result["log_tail"]

    project = load_project(sim.param_file, workdir=tmp_path / "depsi_work")
    selection = run_selection(project)
    n_ifgs = project.Nifgs

    psc_oct, ph_oct, h2_oct = read_points(workdir / "bench_psc_2orig_sel1.raw", n_ifgs, psp=False)
    psp_oct, _, _ = read_points(workdir / "bench_psp_2orig_sel1.raw", n_ifgs, psp=True)
    psc_dep, ph_dep, h2_dep = read_points(project.path("_psc_2orig_sel1.raw"), n_ifgs, psp=False)
    psp_dep, _, _ = read_points(project.path("_psp_2orig_sel1.raw"), n_ifgs, psp=True)

    assert int(selection.Npsc[0]) == result["Npsc"]
    assert int(selection.Npsp[0]) == result["Npsp"]
    psc = compare_points(psc_dep, ph_dep, h2_dep, psc_oct, ph_oct, h2_oct)
    assert psc["exact"], psc
    assert psc["phase_max_abs_rad"] < 1e-5  # float32 raster storage
    assert psc["h2ph_max_abs"] == 0.0

    psp = compare_points(
        psp_dep, np.zeros((psp_dep.shape[0], n_ifgs)), np.zeros((psp_dep.shape[0], n_ifgs)),
        psp_oct, np.zeros((psp_oct.shape[0], n_ifgs)), np.zeros((psp_oct.shape[0], n_ifgs)),
    )
    assert psp["exact"], psp
