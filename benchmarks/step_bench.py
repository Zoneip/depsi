"""Step-3 microbenchmark: selection -> network -> VCE -> estimation -> unwrapping.

This script is the W0 baseline for the vectorisation work.  It uses the same
synthetic project as ``benchmarks/make_scenarios.py`` (``<scenario>/sim``) and
reports wall time per stage plus the process peak RSS.  With ``--profile`` it
also records the top cumulative functions from ``cProfile``.

Example::

    PYTHONPATH=src .venv/bin/python benchmarks/make_scenarios.py --out /tmp/dpbench
    PYTHONPATH=src .venv/bin/python benchmarks/step_bench.py \
        --scenario /tmp/dpbench --json benchmarks/results/w0-step3.json
"""

from __future__ import annotations

import argparse
import cProfile
import json
import pstats
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))

from depsipy.core.estimation import estimate_arcs  # noqa: E402
from depsipy.io.params import load_params  # noqa: E402
from depsipy.io.project import read_input_file, sensor_constants  # noqa: E402
from depsipy.steps.network import (  # noqa: E402
    estimate_variance_components,
    form_network,
)
from depsipy.steps.selection import select_psc_psp  # noqa: E402
from depsipy.steps.unwrapping import unwrap_network  # noqa: E402


def _selection_kwargs(sim: Path):
    prm = load_params(sim / "param_file.txt")
    st = read_input_file(prm.str("input_file"))
    kw = dict(
        slc_files=st.slc_paths(),
        ifg_files=[st.resolve(f) for f in st.ifg],
        h2ph_files=[st.resolve(f) for f in st.h2ph],
        Nlines=st.crop_final[1] - st.crop_final[0] + 1,
        Npixels=st.crop_final[3] - st.crop_final[2] + 1,
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
    return prm, st, kw


def _psc_arrays(result, z: int = 0):
    cand = result.psc[z]
    if not cand:
        raise RuntimeError("selection produced no PSC; step-3 benchmark needs at least one")
    az = np.array([c.az for c in cand], dtype=float)
    r = np.array([c.r for c in cand], dtype=float)
    phase = np.stack([c.phase for c in cand])
    h2ph = np.stack([c.h2ph for c in cand])
    return az, r, phase, h2ph


def run_step3(scenario: Path) -> dict:
    sim = scenario / "sim"
    prm, st, sel_kw = _selection_kwargs(sim)
    m2ph, sat_vel, _ = sensor_constants(prm.str("sensor"))
    Nifgs = st.Nifgs
    method = prm.str("ps_method", "boot")
    model = int(prm.ints("psc_model", [1])[0])
    sig2_est = None
    models = prm.ints("psc_model", [1])
    std_param = prm.floats("std_param")

    timings = {}
    t0 = time.perf_counter()
    selection = select_psc_psp(**sel_kw)
    timings["selection"] = time.perf_counter() - t0

    az, r, phase, h2ph = _psc_arrays(selection)
    t0 = time.perf_counter()
    network = form_network(
        az,
        r,
        phase,
        h2ph,
        az_spacing=sel_kw["az_spacing"],
        r_spacing=sel_kw["r_spacing"],
        max_arc_length=prm.float("max_arc_length", 1000.0),
        method=prm.str("network_method", "delaunay"),
        Ncon=prm.int("Ncon", 8),
    )
    timings["network"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    try:
        sig2_est = estimate_variance_components(
            network.arcs,
            network.dphase,
            network.dh2ph,
            Nifgs=Nifgs,
            Btemp=st.Btemp,
            Bdop=st.Bdop,
            std_param=std_param,
            model=model,
            m2ph=m2ph,
            sat_vel=sat_vel,
            method=method,
        )
        timings["vce"] = time.perf_counter() - t0
    except ValueError:
        # A tiny synthetic network may not have enough independent arcs.
        timings["vce"] = float("nan")

    t0 = time.perf_counter()
    est = estimate_arcs(
        network.dphase,
        Nifgs,
        network.dh2ph,
        st.Btemp,
        st.Bdop,
        std_param,
        sig2_est,
        models,
        m2ph=m2ph,
        sat_vel=sat_vel,
        method=method,
    )
    timings["estimate"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    unwrap_network(
        network.arcs,
        network.dphase,
        est.acheck,
        est.varfac,
        network.Nnodes,
        ps_az=az,
        ps_r=r,
        Nlines=sel_kw["Nlines"],
        Npixels=sel_kw["Npixels"],
        az_spacing=sel_kw["az_spacing"],
        r_spacing=sel_kw["r_spacing"],
    )
    timings["unwrap"] = time.perf_counter() - t0

    return dict(
        schema=1,
        scenario=str(scenario),
        Nlines=sel_kw["Nlines"],
        Npixels=sel_kw["Npixels"],
        Nslc=st.Nslc,
        Nifgs=Nifgs,
        Npsc=len(az),
        Narcs=network.arcs.shape[0],
        stages={k: round(v, 6) for k, v in timings.items()},
        peak_mem_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1),
    )


def _profile(scenario: Path) -> list[dict]:
    prof = cProfile.Profile()
    prof.enable()
    run_step3(scenario)
    prof.disable()
    stats = pstats.Stats(prof)
    out = []
    for func, (_cc, _nc, _tt, cumtime, _callers) in stats.stats.items():
        out.append(dict(function=f"{func[0]}:{func[1]}:{func[2]}", cumtime=round(cumtime, 4)))
    out.sort(key=lambda item: item["cumtime"], reverse=True)
    return out[:20]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scenario", default="/tmp/dpbench")
    ap.add_argument("--json", default=None)
    ap.add_argument("--no-profile", action="store_true")
    args = ap.parse_args(argv)

    scenario = Path(args.scenario)
    if not (scenario / "sim" / "param_file.txt").is_file():
        raise SystemExit(f"{scenario}/sim/param_file.txt not found; run make_scenarios.py first")

    report = run_step3(scenario)
    if not args.no_profile:
        report["profile"] = _profile(scenario)
    text = json.dumps(report, ensure_ascii=False, indent=2)
    if args.json:
        path = Path(args.json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
