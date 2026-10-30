"""Same-data parity run: Octave (MATLAB sources) vs depsipy on one synthetic project.

The MATLAB toolbox is not available on this machine; Octave 9.4 + the shims in
``benchmarks/octave/shim`` execute the *unmodified* ``ps_selection.m`` instead.  This
script generates (or reuses) one DePSI-format project, runs both implementations
on exactly the same rasters/parameters, and compares

* the number of PSC / PSP candidates,
* the PSC and PSP ``(az, r)`` lists (as sets and as sorted lists),
* the PSC phases / h2ph factors at the shared points.

Usage::

    PYTHONPATH=src .venv/bin/python scripts/octave_parity.py \
        --toolbox workspace/ps_toolbox_v1.7.8.8_hhf --out /tmp/dp_parity

    # or reuse an existing project:
    ... --project /path/to/param_file.txt

Exit code 0 means the comparison ran; ``"match": true/false`` says whether the
PSC/PSP sets agree.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))

from depsipy.io.params import load_params  # noqa: E402
from depsipy.io.project import read_input_file  # noqa: E402
from depsipy.simulate import generate_project  # noqa: E402


def write_selection_input(param_file, out_dir: Path) -> Path:
    """``sel_input.mat`` for ``ps_selection.m`` built from the project's own files."""
    import scipy.io as sio

    prm = load_params(param_file)
    stack = read_input_file(prm.str("input_file"))
    path = out_dir / "sel_input.mat"
    sio.savemat(
        path,
        {
            "slc": np.array([stack.resolve(f) for f in stack.slc], dtype=object).reshape(-1, 1),
            "ifg": np.array([stack.resolve(f) for f in stack.ifg], dtype=object).reshape(-1, 1),
            "h2ph": np.array([stack.resolve(f) for f in stack.h2ph], dtype=object).reshape(-1, 1),
            "Btemp": stack.Btemp,
            "crop_in": stack.crop_in.astype(float),
            "crop": np.array(stack.crop_final, dtype=float),
            "slc_selection": np.arange(1, stack.Nslc + 1, dtype=float),
            "calfactors": np.ones(stack.Nslc),
            "psc_selection_method": prm.str("psc_selection_method"),
            "psc_selection_gridsize": prm.float("psc_selection_gridsize"),
            "psc_threshold": prm.float("psc_threshold"),
            "Npsc_selections": prm.int("Npsc_selections"),
            "psp_selection_method": prm.str("psp_selection_method"),
            "psp_threshold1": prm.float("psp_threshold1"),
            "gamma_threshold": prm.float("gamma_threshold"),
            "psc_distribution": prm.str("psc_distribution"),
            "livetime_threshold": prm.float("livetime_threshold"),
            "peak_tolerance": prm.float("peak_tolerance"),
            "max_mem_buffer": 50e7,
            "az_spacing": prm.float("az_spacing"),
            "r_spacing": prm.float("r_spacing"),
            "Nlines": float(stack.crop_final[1] - stack.crop_final[0] + 1),
            "Npixels": float(stack.crop_final[3] - stack.crop_final[2] + 1),
        },
    )
    return path


def run_octave_selection(toolbox: Path, scenario_dir: Path, workdir: Path, *, octave="octave") -> dict:
    """Run ``bench_selection.m`` and return the parsed ``SEL_*`` lines plus timings."""
    script = REPO / "benchmarks/octave/bench_selection.m"
    proc = subprocess.run(
        [octave, "--no-gui", "--quiet", str(script), str(toolbox), str(scenario_dir), str(workdir)],
        capture_output=True, text=True, timeout=3600,
    )
    out = proc.stdout + proc.stderr
    result = {"returncode": proc.returncode, "log_tail": out[-2000:]}
    for line in out.splitlines():
        if line.startswith("SEL_NPSC "):
            result["Npsc"] = int(float(line.split()[1]))
        elif line.startswith("SEL_NPSP "):
            result["Npsp"] = int(float(line.split()[1]))
        elif line.startswith("SEL_TIME "):
            result["time_s"] = float(line.split()[1])
    return result


def read_points(path: Path, n_ifgs: int, *, psp: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(az_r, phase, h2ph)`` of ``_psc_2orig_sel*.raw`` / ``_psp_2orig_sel*.raw``."""
    empty = np.zeros((0, 2), dtype=np.int64)
    if not path.exists():
        return empty, np.zeros((0, n_ifgs)), np.zeros((0, n_ifgs))
    raw = np.fromfile(path, dtype=np.float64)
    if raw.size == 0:
        return empty, np.zeros((0, n_ifgs)), np.zeros((0, n_ifgs))
    cols = (2 * n_ifgs + 4) if psp else (2 * n_ifgs + 7)
    if raw.size % cols:
        raise ValueError(f"{path}: {raw.size} doubles is not a multiple of {cols} columns")
    data = raw.reshape(-1, cols)
    # 0-based columns: PSP = [cell_az cell_r az r phase h2ph]; PSC = [id group grid_az grid_r az r phase h2ph ...]
    az_col, r_col = (2, 3) if psp else (4, 5)
    phase = data[:, r_col + 1 : r_col + 1 + n_ifgs]
    h2ph = data[:, r_col + 1 + n_ifgs : r_col + 1 + 2 * n_ifgs]
    return data[:, [az_col, r_col]].astype(np.int64), phase, h2ph


def compare_points(points_a, phase_a, h2ph_a, points_b, phase_b, h2ph_b) -> dict:
    """Set comparison of the ``(az, r)`` lists plus phase/h2ph agreement."""
    sa = {tuple(row): i for i, row in enumerate(points_a.tolist())}
    sb = {tuple(row): i for i, row in enumerate(points_b.tolist())}
    common = sorted(set(sa) & set(sb))
    out = {
        "count_a": len(sa),
        "count_b": len(sb),
        "common": len(common),
        "only_a": len(set(sa) - set(sb)),
        "only_b": len(set(sb) - set(sa)),
        "exact": set(sa) == set(sb),
    }
    if common:
        ia = np.array([sa[k] for k in common])
        ib = np.array([sb[k] for k in common])
        dphase = (phase_a[ia] - phase_b[ib] + np.pi) % (2 * np.pi) - np.pi
        out["phase_max_abs_rad"] = float(np.max(np.abs(dphase)))
        out["h2ph_max_abs"] = float(np.max(np.abs(h2ph_a[ia] - h2ph_b[ib])))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--toolbox", default="workspace/ps_toolbox_v1.7.8.8_hhf")
    ap.add_argument("--out", default="/tmp/dp_parity")
    ap.add_argument("--project", default=None, help="param_file.txt of an existing project")
    ap.add_argument("--lines", type=int, default=60)
    ap.add_argument("--pixels", type=int, default=60)
    ap.add_argument("--nslc", type=int, default=15)
    ap.add_argument("--seed", type=int, default=2)
    ap.add_argument("--ps_fraction", type=float, default=0.15)
    ap.add_argument("--atmo", type=float, default=0.2)
    ap.add_argument("--defo", type=float, default=0.03)
    ap.add_argument("--gridsize", type=float, default=5.0)
    ap.add_argument("--octave", default="octave")
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)

    octave = shutil.which(args.octave)
    if octave is None:
        print(f"error: {args.octave!r} not found on PATH", file=sys.stderr)
        return 2
    toolbox = (REPO / args.toolbox).resolve()
    if not (toolbox / "main/ps_selection.m").exists():
        print(f"error: MATLAB toolbox not found at {toolbox}", file=sys.stderr)
        return 2

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.project:
        param_file = Path(args.project).resolve()
    else:
        proj = generate_project(
            out / "sim",
            Nlines=args.lines, Npixels=args.pixels, Nslc=args.nslc, seed=args.seed,
            ps_fraction=args.ps_fraction, atmo_amp=args.atmo, defo_rate=args.defo,
            gridsize=args.gridsize,
        )
        param_file = proj.param_file

    scenario = out / "octave_input"
    scenario.mkdir(parents=True, exist_ok=True)
    write_selection_input(param_file, scenario)
    workdir = out / "octave_work"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    # depsipy on the same project
    from depsipy.pipeline import load_project, run_selection

    project = load_project(param_file, workdir=out / "depsi_work")
    selection = run_selection(project)
    n_ifgs = project.Nifgs

    octave_result = run_octave_selection(toolbox, scenario, workdir, octave=octave)
    if octave_result["returncode"] != 0:
        print(json.dumps(octave_result, indent=2))
        return 1

    psc_oct, ph_oct, h2_oct = read_points(workdir / "bench_psc_2orig_sel1.raw", n_ifgs, psp=False)
    psp_oct, _, _ = read_points(workdir / "bench_psp_2orig_sel1.raw", n_ifgs, psp=True)
    psc_dep, ph_dep, h2_dep = read_points(project.path("_psc_2orig_sel1.raw"), n_ifgs, psp=False)
    psp_dep, _, _ = read_points(project.path("_psp_2orig_sel1.raw"), n_ifgs, psp=True)

    summary = {
        "project": str(param_file),
        "Nslc": project.Nslc,
        "Nifgs": n_ifgs,
        "octave": {k: octave_result.get(k) for k in ("Npsc", "Npsp", "time_s")},
        "depsi": {"Npsc": int(selection.Npsc[0]), "Npsp": int(selection.Npsp[0])},
        "psc": compare_points(psc_dep, ph_dep, h2_dep, psc_oct, ph_oct, h2_oct),
        "psp": compare_points(psp_dep, np.zeros((psp_dep.shape[0], n_ifgs)),
                              np.zeros((psp_dep.shape[0], n_ifgs)),
                              psp_oct, np.zeros((psp_oct.shape[0], n_ifgs)),
                              np.zeros((psp_oct.shape[0], n_ifgs))),
        "match": (
            int(selection.Npsc[0]) == octave_result.get("Npsc")
            and int(selection.Npsp[0]) == octave_result.get("Npsp")
        ),
    }
    text = json.dumps(summary, indent=2)
    print(text)
    if args.json:
        Path(args.json).write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
