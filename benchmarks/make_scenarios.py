"""Build the benchmark scenarios shared by Octave (MATLAB code) and depsipy.

Two scenarios are produced under ``<out>/``:

``sim/``            a full DePSI-format project (raw rasters + input_file + param_file)
                    used to time the selection step end to end.
``est_input.mat``   the arrays needed by ``ps_ils_and_bootstrap_estimation.m``
``est_input.npz``   the same arrays for depsipy

The estimation arrays are *identical* for both sides, so the timing difference is
implementation, not data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from depsipy.simulate import generate_project


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/dpbench")
    ap.add_argument("--lines", type=int, default=400)
    ap.add_argument("--pixels", type=int, default=400)
    ap.add_argument("--nslc", type=int, default=20)
    ap.add_argument("--arcs", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # ---- scenario 1: selection -------------------------------------------
    proj = generate_project(
        out / "sim",
        Nlines=args.lines,
        Npixels=args.pixels,
        Nslc=args.nslc,
        seed=args.seed,
        ps_fraction=0.10,
    )

    # ---- scenario 2: per-arc estimation ----------------------------------
    rng = np.random.default_rng(args.seed + 1)
    Nifgs = args.nslc - 1
    Btemp = np.sort(rng.uniform(-1.0, 1.0, Nifgs))
    Bdop = np.round(rng.normal(0, 10, Nifgs), 1)
    bperp = np.linspace(-600, 600, Nifgs)
    h2ph = -2 * bperp / (0.055517120823757024 * 800e3 * np.sin(np.deg2rad(23)))
    # per-arc mean height-to-phase factor with a little scatter
    h2ph_arcs = np.tile(h2ph, (args.arcs, 1)) * rng.uniform(0.9, 1.1, (args.arcs, 1))
    dphase = np.mod(rng.normal(0, 0.7, (args.arcs, Nifgs)) + np.pi, 2 * np.pi) - np.pi
    std_param = np.array([30.0, 0.005, 1.0, 0.02, 0.01, 0.01, 0.005])
    psc_model = np.array([1.0])
    final_model = np.array([2.0])

    import scipy.io as sio

    sio.savemat(
        out / "est_input.mat",
        {
            "h2ph": h2ph_arcs,
            "Btemp": Btemp,
            "Bdop": Bdop,
            "dphase": dphase,
            "model": psc_model,
            "std_param": std_param,
            "sig2_est": np.zeros((0, 0)),
            "breakpoint": np.zeros((0, 0)),
            "breakpoint2": np.zeros((0, 0)),
        },
    )
    np.savez(
        out / "est_input.npz",
        h2ph=h2ph_arcs,
        Btemp=Btemp,
        Bdop=Bdop,
        dphase=dphase,
        std_param=std_param,
    )

    # ---- scenario 1b: the selection arrays, prepared for MATLAB -----------
    from depsipy.io.params import load_params
    from depsipy.io.project import read_input_file

    prm = load_params(proj.param_file)
    st = read_input_file(prm.str("input_file"))
    sio.savemat(
        out / "sel_input.mat",
        {
            # column cells: MATLAB's ps_selection uses size(filenames_slc,1) as
            # the number of SLCs, so a row vector would be read as a single scene
            "slc": np.array([st.resolve(f) for f in st.slc], dtype=object).reshape(-1, 1),
            "ifg": np.array([st.resolve(f) for f in st.ifg], dtype=object).reshape(-1, 1),
            "h2ph": np.array([st.resolve(f) for f in st.h2ph], dtype=object).reshape(-1, 1),
            "Btemp": st.Btemp,
            "crop_in": st.crop_in.astype(float),
            "crop": np.array(st.crop_final, dtype=float),
            "slc_selection": np.arange(1, st.Nslc + 1, dtype=float),
            "calfactors": np.ones(st.Nslc),
            "psc_selection_method": prm.str("psc_selection_method"),
            "psc_selection_gridsize": prm.float("psc_selection_gridsize"),
            "psc_threshold": prm.float("psc_threshold"),
            "Npsc_selections": prm.int("Npsc_selections"),
            "psp_selection_method": prm.str("psp_selection_method"),
            "psp_threshold1": prm.float("psp_threshold1"),
            "gamma_threshold": prm.float("gamma_threshold"),
            "psc_distribution": prm.str("psc_distribution"),
            "livetime_threshold": prm.float("livetime_threshold"),
            "peak_tolerance": 0.7,
            "max_mem_buffer": 50e7,
            "az_spacing": prm.float("az_spacing"),
            "r_spacing": prm.float("r_spacing"),
            "Nlines": float(args.lines),
            "Npixels": float(args.pixels),
        },
    )

    meta = dict(
        lines=args.lines,
        pixels=args.pixels,
        nslc=args.nslc,
        nifgs=Nifgs,
        arcs=args.arcs,
        seed=args.seed,
        sim=dict(workdir=str(proj.workdir), param_file=str(proj.param_file)),
    )
    (out / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
