"""Round-trip the numeric MATLAB interchange without a MATLAB dependency."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.io import savemat

SCRIPT = Path(__file__).resolve().parents[1]/'scripts/readiness.py'


def test_mat_to_npz_conversion(tmp_path):
    source, target = tmp_path/'source.mat', tmp_path/'target.npz'
    savemat(source, {'ids': np.array([[2,3],[4,5]]), 'epochs': np.array([[20250101,20250102]]), 'v': [[1.,2.],[3.,4.]]})
    command = [sys.executable, str(SCRIPT), 'from-mat', str(source), str(target), '--field', 'v']
    p = subprocess.run(command, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    with np.load(target) as data:
        np.testing.assert_array_equal(data['ids'], [[2,3],[4,5]])
        np.testing.assert_array_equal(data['epochs'], [20250101,20250102])
        np.testing.assert_array_equal(data['v'], [[1.,2.],[3.,4.]])
    assert subprocess.run(command, capture_output=True).returncode == 2


def test_nan_only_is_not_a_valid_comparison(tmp_path):
    spec = importlib.util.spec_from_file_location('readiness_nan', SCRIPT)
    tool = importlib.util.module_from_spec(spec); spec.loader.exec_module(tool)
    a, b = tmp_path/'a.npz', tmp_path/'b.npz'
    for p in (a,b):
        np.savez(p, ids=[1], v=[np.nan])
    meta = tmp_path/'meta.json'
    meta.write_text(json.dumps({'units':{'v':'m'}, 'reference':'same', 'phase_direction':'slave-minus-master'}))
    report = tool.compare(a,b,meta,meta,['v'])
    assert report['status'] == 'blocked'
