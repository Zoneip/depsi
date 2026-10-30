"""Small, filesystem-local checks for preparation tooling, not cloud execution."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'readiness.py'
spec = importlib.util.spec_from_file_location('readiness', SCRIPT)
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)


@pytest.fixture
def project(tmp_path):
    for name in ('slave.raw', 'master.raw', 'ifg.raw'):
        np.ones((2, 3), dtype=np.complex64).tofile(tmp_path/name)
    np.ones((2, 3), dtype=np.float32).tofile(tmp_path/'h2ph.raw')
    text = (f"2 20250102 {tmp_path/'slave.raw'} {tmp_path/'ifg.raw'} {tmp_path/'h2ph.raw'} /original/output 10 0.1 0 1 2 1 3\n"
            f"1 20250101 {tmp_path/'master.raw'} dummy dummy dummy 0 0 0 1 2 1 3\n")
    (tmp_path/'input.txt').write_text(text)
    params = tmp_path/'params.txt'
    params.write_text(f"project_id = 'tiny'\ninput_file = '{tmp_path/'input.txt'}'\n"
                      "sensor = 'spacety'\norbit = 'asc'\naz_spacing = 1.8\nr_spacing = 0.62\n"
                      "processing_groups = [1:8]\ndefo_method = []\nps_method = 'boot'\n")
    return params


def codes(report):
    return {i['code'] for i in report['issues']}


def test_preflight_readonly(project):
    before = {p: p.read_bytes() for p in project.parent.iterdir()}
    report = tool.preflight(project, sample_pixels=3)
    assert report['status'] == 'pass'
    assert report['row_count'] == 2
    assert 'group6_skipped' in codes(report)
    assert all(p.read_bytes() == data for p, data in before.items())


def test_missing_and_bad_size(project):
    (project.parent/'slave.raw').unlink()
    (project.parent/'h2ph.raw').write_bytes(b'bad')
    report = tool.preflight(project)
    assert report['status'] == 'blocked'
    assert {'raster_missing', 'raster_size'} <= codes(report)


def test_duplicate_master(project):
    path = project.parent/'input.txt'
    text = path.read_text()
    path.write_text(text + text.splitlines()[1] + '\n')
    report = tool.preflight(project)
    assert {'master_count', 'duplicate_acquisition'} <= codes(report)


def test_step9_is_pending_not_silent_success(project):
    project.write_text(project.read_text().replace('[1:8]', '[1:9]'))
    assert 'step_mapping_pending' in codes(tool.preflight(project))
    assert tool.preflight(project, baseline_policy=True)['status'] == 'blocked'


def test_path_prefix_components(tmp_path):
    assert tool.resolve_path('/old/x', tmp_path, {'/old': str(tmp_path)}) == tmp_path/'x'
    assert str(tool.resolve_path('/older/x', tmp_path, {'/old': '/new'})) == '/older/x'


def test_sampling_nan_is_warning(project):
    a = np.ones(6, dtype=np.float32); a[0] = np.nan
    a.tofile(project.parent/'h2ph.raw')
    assert 'nonfinite_sample' in codes(tool.preflight(project, sample_pixels=6))


@pytest.fixture
def pair(tmp_path):
    meta = dict(units={'v': 'm'}, reference={'az': 1, 'r': 1}, phase_direction='slave-minus-master')
    lm, rm = tmp_path/'left.json', tmp_path/'right.json'
    lm.write_text(json.dumps(meta)); rm.write_text(json.dumps(meta))
    left, right = tmp_path/'a.npz', tmp_path/'b.npz'
    np.savez(left, ids=[1, 2], epochs=['20250101', '20250102'], v=[[1., 2.], [3., np.nan]])
    np.savez(right, ids=[2, 1], epochs=['20250102', '20250101'], v=[[np.nan, 3.], [2., 1.]])
    return left, right, lm, rm


def test_compare_aligns_points_and_dates(pair):
    r = tool.compare(*pair, ['v'])
    assert r['status'] == 'pass'
    assert r['fields']['v']['max_abs'] == 0


def test_compare_identity_set_difference(pair):
    np.savez(pair[1], ids=[2, 3], epochs=['20250101', '20250102'], v=[[3., np.nan], [8., 9.]])
    r = tool.compare(*pair, ['v'])
    assert r['status'] == 'mismatch'
    assert r['identity_sets']['left_only'] == 1


def test_compare_reference_difference_blocked(pair):
    d = json.loads(pair[3].read_text()); d['reference'] = {'az': 2, 'r': 2}
    pair[3].write_text(json.dumps(d))
    assert tool.compare(*pair, ['v'])['status'] == 'blocked'


def test_compare_nan_and_inf_sign(pair):
    np.savez(pair[0], ids=[1], v=[np.inf])
    np.savez(pair[1], ids=[1], v=[-np.inf])
    assert tool.compare(*pair, ['v'])['status'] == 'mismatch'
    np.savez(pair[1], ids=[1], v=[np.nan])
    assert tool.compare(*pair, ['v'])['status'] == 'mismatch'


def test_duplicate_ids_rejected(pair):
    np.savez(pair[1], ids=[1, 1], v=[1., 2.])
    with pytest.raises(ValueError, match='duplicate'):
        tool.compare(*pair, ['v'])


def test_cli_failure_is_nonzero_and_json(tmp_path):
    p = subprocess.run([sys.executable, str(SCRIPT), 'preflight', str(tmp_path/'missing')], capture_output=True, text=True)
    assert p.returncode == 2
    assert json.loads(p.stdout)['status'] == 'blocked'


def mat_args(project):
    entry = project.parent/'depsi.m'
    entry.write_text('function depsi(params)\nend\n')
    return SimpleNamespace(entry=str(entry), params=str(project), run_dir=str(project.parent/'run'),
                           matlab='matlab', execute=False, profile=True, timeout=10)


def test_matlab_plan_is_not_execution_and_isolates_outputs(project):
    args = mat_args(project)
    original = project.read_bytes()
    report = tool.matlab_plan(args)
    assert report['status'] == 'prepared' and report['executed'] is False
    assert report['command'][1] == '-batch'
    cfg = report['config']
    assert 'shim' not in report['command'][2]
    assert Path(cfg['param_file']).read_text() != project.read_text()
    staged_input = Path(args.run_dir)/'input.staged.txt'
    assert '/original/output' not in staged_input.read_text()
    assert project.read_bytes() == original
    with pytest.raises(FileExistsError):
        tool.matlab_plan(args)


def test_missing_matlab_propagates_failure(project):
    args = mat_args(project); args.execute = True
    args.matlab = '/no/such/matlab'
    report = tool.matlab_plan(args)
    assert report['status'] == 'failed'
    assert report['executed'] is True and report['returncode'] != 0


def test_zero_exit_without_result_is_not_success(project, monkeypatch):
    args = mat_args(project); args.execute = True
    monkeypatch.setattr(tool.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0))
    assert tool.matlab_plan(args)['status'] == 'failed'


def test_matlab_plan_resolves_relative_rasters(project):
    """Relative names are pinned to the input-file directory and recorded."""
    text = (project.parent/'input.txt').read_text()
    (project.parent/'input.txt').write_text('\n'.join(
        ' '.join([c[0], c[1], Path(c[2]).name, Path(c[3]).name if c[3] != 'dummy' else 'dummy',
                  Path(c[4]).name if c[4] != 'dummy' else 'dummy', c[5]] + c[6:])
        for c in (line.split() for line in text.splitlines())
    ) + '\n')
    report = tool.matlab_plan(mat_args(project))
    assert report['status'] == 'prepared'
    resolved = report['config']['resolved_relative_paths']
    assert resolved and all(Path(r['staged']).is_file() for r in resolved)
    staged = Path(report['config']['param_file']).parent/'input.staged.txt'
    assert str(project.parent) in staged.read_text()


def test_matlab_plan_refuses_missing_raster(project):
    (project.parent/'slave.raw').unlink()
    args = mat_args(project)
    text = (project.parent/'input.txt').read_text().replace(
        str(project.parent/'slave.raw'), 'slave.raw')
    (project.parent/'input.txt').write_text(text)
    report = tool.matlab_plan(args)
    assert report['status'] == 'blocked'


def test_compare_cli_exit_codes(pair, tmp_path):
    out = tmp_path/'cmp.json'
    ok = subprocess.run([sys.executable, str(SCRIPT), 'compare', str(pair[0]), str(pair[1]),
                         '--meta-left', str(pair[2]), '--meta-right', str(pair[3]),
                         '--field', 'v', '--output', str(out)], capture_output=True, text=True)
    assert ok.returncode == 0, ok.stdout + ok.stderr
    np.savez(pair[1], ids=[2, 1], epochs=['20250102', '20250101'], v=[[np.nan, 3.5], [2., 1.]])
    bad = subprocess.run([sys.executable, str(SCRIPT), 'compare', str(pair[0]), str(pair[1]),
                          '--meta-left', str(pair[2]), '--meta-right', str(pair[3]),
                          '--field', 'v'], capture_output=True, text=True)
    assert bad.returncode == 1
    assert json.loads(bad.stdout)['status'] == 'mismatch'
