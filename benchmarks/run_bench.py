"""Run the depsipy side of the benchmarks and (optionally) compare with Octave.

    python benchmarks/run_bench.py --scenario /tmp/dpbench
    python benchmarks/run_bench.py --scenario /tmp/dpbench --with-octave workspace/ps_toolbox_v1.7.8.8_hhf

The Octave side is driven from ``benchmarks/octave/bench_*.m``; when
``--with-octave`` is given this script also times the MATLAB implementation and
prints a side-by-side table.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))

from depsipy.core.estimation import estimate_arcs  # noqa: E402
from depsipy.io.params import load_params  # noqa: E402
from depsipy.io.project import read_input_file  # noqa: E402
from depsipy.steps.selection import select_psc_psp  # noqa: E402

M2PH = -4 * np.pi / 0.055517120823757024
SAT_VEL = 7690.0


def bench_estimation(scen: Path, reps: int = 5) -> dict:
    d = np.load(scen / "est_input.npz")
    h2ph, Btemp, Bdop, dphase, std = d["h2ph"], d["Btemp"], d["Bdop"], d["dphase"], d["std_param"]
    n = dphase.shape[1]
    kw = dict(m2ph=M2PH, sat_vel=SAT_VEL, method="boot")
    estimate_arcs(dphase[:20], n, h2ph[:20], Btemp, Bdop, std, None, [1], **kw)  # JIT warm-up

    times = []
    for _ in range(reps):
        t = time.perf_counter()
        res = estimate_arcs(dphase, n, h2ph, Btemp, Bdop, std, None, [1], **kw)
        times.append(time.perf_counter() - t)

    out = dict(
        arcs=int(dphase.shape[0]),
        nifgs=int(n),
        time_s=round(min(times), 4),
        reps=[round(t, 4) for t in times],
        throughput_arcs_s=int(dphase.shape[0] / min(times)),
        median_varfac=float(np.median(res.varfac[np.isfinite(res.varfac)])),
        median_ens_coh=float(np.median(res.ens_coh[np.isfinite(res.ens_coh)])),
    )
    ref = scen / "ref3000.csv"
    if ref.exists():
        r = np.loadtxt(ref, delimiter=",")
        out["session"] = dict(
            acheck_maxdiff=float(np.abs(res.acheck - r[:, 2:]).max()),
            varfac_maxdiff=float(np.abs(res.varfac - r[:, 0]).max()),
            ens_coh_maxdiff=float(np.abs(res.ens_coh - r[:, 1]).max()),
        )
    return out


def bench_selection(scen: Path, reps: int = 3) -> dict:
    sim = scen / "sim"
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
    select_psc_psp(**kw)  # page-cache warm-up
    times = []
    for _ in range(reps):
        t = time.perf_counter()
        res = select_psc_psp(**kw)
        times.append(time.perf_counter() - t)
    return dict(
        lines=kw["Nlines"],
        pixels=kw["Npixels"],
        nslc=st.Nslc,
        time_s=round(min(times), 4),
        npsc=int(res.Npsc[0]),
        npsp=int(res.Npsp[0]),
    )


def run_octave(script: str, toolbox: Path, scen: Path, extra: Path | None = None) -> dict:
    cmd = [
        "octave", "--no-gui", "--quiet", str(HERE / "octave" / script),
        str(toolbox.resolve()), str(scen.resolve()),
    ]
    if extra is not None:
        cmd.append(str(extra))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600,
                          cwd=str(HERE / "octave"))
    out = {}
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isupper() and parts[0].startswith(("EST_", "SEL_", "PEAK_")):
            key, val = parts
            try:
                out[key.lower()] = float(val)
            except ValueError:
                pass
    out["_returncode"] = proc.returncode
    if proc.returncode != 0:
        out["_stderr_tail"] = proc.stderr.strip().splitlines()[-3:]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="/tmp/dpbench")
    ap.add_argument("--with-octave", default=None,
                    help="toolbox directory of the MATLAB version to benchmark")
    ap.add_argument("--json", default=None, help="write the results here")
    args = ap.parse_args(argv)

    scen = Path(args.scenario)
    results = {"depsipy": {"estimation": bench_estimation(scen),
                          "selection": bench_selection(scen)}}

    if args.with_octave:
        toolbox = Path(args.with_octave)
        tmp = Path(tempfile.mkdtemp(prefix="dpsel_"))
        try:
            m_est = run_octave("bench_estimation.m", toolbox, scen)
            m_sel = run_octave("bench_selection.m", toolbox, scen, tmp)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        results["matlab_via_octave"] = {"toolbox": str(toolbox), "estimation": m_est, "selection": m_sel}
        py_est = results["depsipy"]["estimation"]["time_s"]
        py_sel = results["depsipy"]["selection"]["time_s"]
        if m_est.get("est_time"):
            results["speedup"] = {
                "estimation": round(m_est["est_time"] / py_est, 1),
                "selection": round(m_sel["sel_time"] / py_sel, 1) if m_sel.get("sel_time") else None,
            }

    print(json.dumps(results, indent=2, ensure_ascii=False))
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
