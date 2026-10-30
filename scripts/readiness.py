"""Preparation CLI (no pipeline/network): environment, preflight, compare, matlab.
Exit: 0 pass/prepared/completed, 1 mismatch/failed, 2 blocked.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata as md
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from depsipy.io.params import parse_param_file


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def emit(report, output=None):
    text = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)
    if output:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x', encoding='utf-8') as f:
            f.write(text + '\n')
    print(text)


def environment():
    packages = {}
    for name in ('numpy', 'scipy', 'numba', 'pytest', 'setuptools'):
        try:
            packages[name] = md.version(name)
        except md.PackageNotFoundError:
            packages[name] = None
    check = subprocess.run([sys.executable, '-m', 'pip', 'check'], capture_output=True, text=True)
    return dict(schema=1, status='pass' if check.returncode == 0 else 'blocked',
                python={'version': platform.python_version(), 'executable': sys.executable},
                platform=platform.platform(), cpu_count=os.cpu_count(), packages=packages,
                depsipy_origin=importlib.util.find_spec('depsipy').origin,
                matlab=shutil.which('matlab'), octave=shutil.which('octave'),
                threads={k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS')},
                pip_check={'returncode': check.returncode, 'output': (check.stdout + check.stderr).strip()},
                unverified=['MATLAB license/toolboxes', 'cloud hardware', 'real project'])


def resolve_path(raw, base, mappings):
    for source in sorted(mappings, key=len, reverse=True):
        if raw == source or raw.startswith(source.rstrip('/') + '/'):
            raw = mappings[source].rstrip('/') + raw[len(source):]
            break
    path = Path(raw).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def preflight(param_file, mappings=None, sample_pixels=0, baseline_policy=False):
    mappings = mappings or {}
    param_file = Path(param_file).resolve()
    issues, files = [], []

    def issue(level, code, detail):
        issues.append(dict(level=level, code=code, detail=str(detail)))

    params = parse_param_file(param_file)
    for name in ('project_id', 'input_file', 'sensor', 'orbit', 'az_spacing', 'r_spacing', 'processing_groups'):
        if name not in params or np.size(params[name]) == 0:
            issue('blocked', 'parameter_missing', name)
    for name in ('az_spacing', 'r_spacing', 'max_mem_buffer'):
        if name in params:
            a = np.asarray(params[name], dtype=float)
            if a.size != 1 or not np.isfinite(a).all() or (a <= 0).any():
                issue('blocked', 'parameter_invalid', name)
    groups = np.atleast_1d(params.get('processing_groups', [])).astype(float)
    if not np.isfinite(groups).all() or (groups != np.floor(groups)).any() or (groups < 1).any():
        issue('blocked', 'groups_invalid', groups.tolist())
    if (groups > 8).any():
        issue('blocked' if baseline_policy else 'warning', 'step_mapping_pending',
              'Reference groups are 1..8; production step9 awaits actual entry.')
    if groups.size and 1 not in groups:
        issue('blocked' if baseline_policy else 'unverified', 'checkpoint_required',
              'Starting after group1 needs a complete isolated state; resume launcher not prepared.')
    if 6 in groups and np.size(params.get('defo_method', [])) == 0:
        issue('warning', 'group6_skipped', 'defo_method=[] skips reference group6 modeling/densification branch.')
    if params.get('ps_method', 'boot') not in ('boot', 'ils', 'perio'):
        issue('blocked', 'method_invalid', params['ps_method'])
    if params.get('ps_method') == 'perio':
        issue('warning', 'python_gap', 'This depsipy snapshot has no periodogram estimator.')
    for name in ('psc_model', 'ps_model', 'final_model'):
        if name in params:
            a = np.atleast_1d(params[name]).astype(float)
            if not a.size or not np.isfinite(a).all() or (a != np.floor(a)).any() or ((a < 1) | (a > 27)).any():
                issue('blocked', 'model_invalid', name)
    if 'std_param' in params:
        a = np.atleast_1d(params['std_param']).astype(float)
        if a.size != 7 or not np.isfinite(a).all() or (a <= 0).any():
            issue('blocked', 'std_param_invalid', 'Expected seven positive finite values; units need metadata.')
    for name in ('do_apriori_sidelobe_mask', 'psp_selection_method', 'dens_check'):
        if name in params:
            issue('unverified', 'branch_support', f'{name}={params[name]}: consult gap matrix, not a Python support claim.')
    rows, source = [], None
    if not isinstance(params.get('input_file'), str):
        issue('blocked', 'input_file_invalid', 'Expected a path string.')
    else:
        source = resolve_path(params['input_file'], param_file.parent, mappings)
        if not source.is_file():
            issue('blocked', 'input_missing', source)
        else:
            for number, line in enumerate(source.read_text(encoding='utf-8').splitlines(), 1):
                if not line.strip() or line.lstrip().startswith('%'):
                    continue
                cells = line.split()
                if len(cells) != 13:
                    issue('blocked', 'column_count', f'line {number}: {len(cells)} != 13')
                    continue
                try:
                    nums = np.asarray([float(v) for v in cells[6:]], dtype=float)
                    if not np.isfinite(nums).all() or (nums[3:] != np.floor(nums[3:])).any():
                        raise ValueError('nonfinite baseline or noninteger crop')
                    crop = nums[3:].astype(int)
                    if min(crop) < 1 or crop[1] < crop[0] or crop[3] < crop[2]:
                        raise ValueError('invalid crop')
                except ValueError as exc:
                    issue('blocked', 'row_invalid', f'line {number}: {exc}')
                    continue
                rows.append((cells, nums, crop))
            if not rows:
                issue('blocked', 'no_valid_rows', source)
    masters = [r for r in rows if r[0][3].lower() == 'dummy']
    if rows and len(masters) != 1:
        issue('blocked', 'master_count', len(masters))
    common = None
    if rows:
        ids = [(r[0][0], r[0][1]) for r in rows]
        if len(set(ids)) != len(ids):
            issue('blocked', 'duplicate_acquisition', 'Repeated orbit/date pair (including master).')
        crops = np.array([r[2] for r in rows])
        common = [int(crops[:, 0].max()), int(crops[:, 1].min()), int(crops[:, 2].max()), int(crops[:, 3].min())]
        if common[0] > common[1] or common[2] > common[3]:
            issue('blocked', 'crop_intersection_empty', common)
        issue('unverified', 'baseline_units', 'Confirm Btemp units/sign, master epoch and dates; no automatic unit guessing.')
        issue('unverified', 'physical_units', 'Confirm wavelength, h2ph sign/units and reference height against master_res.')
        bt = [r[1][1] for r in rows if r[0][3].lower() != 'dummy']
        if bt != sorted(bt):
            issue('warning', 'baseline_unsorted', 'Comparison must align acquisitions, not file order.')
    seen = set()
    for cells, nums, crop in rows:
        for col, dtype, size in ((2, 'complex64', 8), (3, 'complex64', 8), (4, 'float32', 4)):
            if cells[col].lower() == 'dummy':
                continue
            path = resolve_path(cells[col], source.parent, mappings)
            expected = int((crop[1]-crop[0]+1)*(crop[3]-crop[2]+1)*size)
            key = (str(path), expected, dtype)
            if key in seen:
                continue
            seen.add(key)
            item = dict(path=str(path), dtype_assumed=dtype, expected_bytes=expected)
            files.append(item)
            if not path.is_file():
                issue('blocked', 'raster_missing', path)
                item['status'] = 'missing'
                continue
            item['bytes'] = path.stat().st_size
            if item['bytes'] != expected:
                issue('blocked', 'raster_size', f'{path}: {item["bytes"]} != {expected}')
                item['status'] = 'size_mismatch'
                continue
            item['status'] = 'size_checked'
            if sample_pixels:
                a = np.memmap(path, dtype=dtype, mode='r')
                idx = np.linspace(0, a.size-1, min(sample_pixels, a.size), dtype=int)
                count = int(np.count_nonzero(~np.isfinite(a[idx])))
                item['sample'] = dict(pixels=len(idx), nonfinite=count)
                if count:
                    issue('warning', 'nonfinite_sample', path)
                del a
    res = params.get('master_res')
    if isinstance(res, str) and res:
        resolved = resolve_path(res, param_file.parent, mappings)
        if not resolved.is_file():
            issue('blocked', 'metadata_missing', resolved)
    else:
        issue('unverified', 'metadata_absent', 'No master_res; geocoding/wavelength not verified.')
    return dict(schema=1, status='blocked' if any(i['level']=='blocked' for i in issues) else 'pass',
                scope='Input structure only, not scientific correctness or pipeline readiness',
                param_file=str(param_file), param_sha256=digest(param_file),
                input_file=str(source) if source else None,
                input_sha256=digest(source) if source and source.is_file() else None,
                row_count=len(rows), common_crop=common, issues=issues, files=files)


def identities(values):
    values = np.asarray(values)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or (values.dtype.kind in 'fc' and not np.isfinite(values).all()):
        raise ValueError('ids must be finite vector or N x D coordinates/endpoints')
    keys = [tuple(row.tolist()) for row in values]
    if len(keys) != len(set(keys)):
        raise ValueError('duplicate identities')
    return keys


def compare(left, right, meta_left, meta_right, fields, atol=1e-8, rtol=1e-7):
    report = dict(schema=1, status='pass', fields={})
    ml, mr = json.loads(Path(meta_left).read_text()), json.loads(Path(meta_right).read_text())
    for name in ('units', 'reference', 'phase_direction'):
        if name not in ml or name not in mr or ml[name] != mr[name]:
            return dict(schema=1, status='blocked', reason=f'metadata missing/different: {name}')
    with np.load(left, allow_pickle=False) as a, np.load(right, allow_pickle=False) as b:
        ia, ib = identities(a['ids']), identities(b['ids'])
        lookup = {k: i for i, k in enumerate(ib)}
        common = [k for k in ia if k in lookup]
        only_a, only_b = set(ia)-set(ib), set(ib)-set(ia)
        report['identity_sets'] = dict(left=len(ia), right=len(ib), common=len(common), left_only=len(only_a), right_only=len(only_b))
        if only_a or only_b or not common:
            report['status'] = 'mismatch'
        aa, bb = [i for i, k in enumerate(ia) if k in lookup], [lookup[k] for k in common]
        has_epochs = 'epochs' in a or 'epochs' in b
        ea = eb = None
        if has_epochs:
            if 'epochs' not in a or 'epochs' not in b:
                return dict(schema=1, status='blocked', reason='epochs missing on one side')
            ka, kb = identities(a['epochs']), identities(b['epochs'])
            if set(ka) != set(kb):
                return dict(schema=1, status='mismatch', reason='epoch sets differ')
            ea, eb = list(range(len(ka))), [kb.index(k) for k in ka]
        report['tolerances'] = dict(atol=atol, rtol=rtol)
        for field in fields:
            if field not in ml['units']:
                raise ValueError(f'unit undeclared: {field}')
            x, y = a[field], b[field]
            if x.shape[0] != len(ia) or y.shape[0] != len(ib):
                raise ValueError(f'{field}: first dimension must match ids')
            x, y = x[aa], y[bb]
            if has_epochs and x.ndim > 1:
                if x.ndim != 2 or y.ndim != 2 or x.shape[1] != len(ea) or y.shape[1] != len(eb):
                    raise ValueError(f'{field}: timeseries must be N x epochs')
                x, y = x[:, ea], y[:, eb]
            if x.shape != y.shape or x.dtype.kind not in 'biufc' or y.dtype.kind not in 'biufc':
                raise ValueError(f'{field}: numeric shapes differ/invalid')
            masks_equal = np.array_equal(np.isnan(x), np.isnan(y)) and np.array_equal(np.isinf(x), np.isinf(y))
            nonfinite_same = bool(
                np.array_equal(np.isnan(x), np.isnan(y))
                and np.array_equal(np.isinf(x), np.isinf(y))
                and np.array_equal(x[np.isinf(x)], y[np.isinf(y)])
            )
            finite = np.isfinite(x) & np.isfinite(y)
            exact = field in ml.get('exact_fields', []) or field in mr.get('exact_fields', []) or (x.dtype.kind in 'biu' and y.dtype.kind in 'biu')
            close = np.equal(x, y) if exact else np.isclose(x, y, atol=atol, rtol=rtol, equal_nan=True)
            close = close | (np.isnan(x) & np.isnan(y))
            ok = masks_equal and nonfinite_same and bool(np.all(close))
            delta = np.abs(x[finite].astype(np.complex128)-y[finite].astype(np.complex128))
            report['fields'][field] = dict(passed=ok, exact=exact, masks_equal=masks_equal,
                nonfinite_equal=nonfinite_same,
                finite_count=int(finite.sum()), nonfinite_count=int((~finite).sum()),
                max_abs=float(delta.max()) if delta.size else None,
                rms=float(np.sqrt(np.mean(delta**2))) if delta.size else None)
            if not masks_equal or not nonfinite_same:
                report['status'] = 'mismatch'
            elif not finite.any():
                report['status'] = 'blocked'
                report['fields'][field]['reason'] = 'No comparable finite values; matching masks alone are not validation.'
            elif not ok and report['status'] != 'blocked':
                report['status'] = 'mismatch'
        report['evidence'] = dict(left_sha256=digest(left), right_sha256=digest(right))
    return report


def matlab_plan(args):
    entry, params = Path(args.entry).resolve(), Path(args.params).resolve()
    if not entry.is_file() or entry.suffix != '.m' or not re.fullmatch(r'[A-Za-z]\w*', entry.stem):
        raise ValueError('entry must be a MATLAB function .m accepting one parameter-file path')
    check = preflight(params, baseline_policy=True)
    if check['status'] != 'pass':
        return dict(schema=1, status='blocked', preflight=check)
    parsed = parse_param_file(params)
    if not Path(parsed['input_file']).is_absolute():
        raise ValueError('MATLAB plan needs absolute input_file; stage a separate parameter copy first')
    input_path = Path(parsed['input_file'])
    for name in ('master_res', 'filename_water_mask', 'slc_selection_input', 'ifg_selection_input'):
        value = parsed.get(name)
        if isinstance(value, str) and value and not Path(value).is_absolute():
            raise ValueError(f'MATLAB plan requires absolute {name}; stage a separate copy first')
    input_lines, relativized = [], []
    for line in input_path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('%'):
            continue
        cells = line.split()
        for i in (2, 3, 4):
            if cells[i].lower() == 'dummy':
                continue
            path_i = Path(cells[i])
            if not path_i.is_absolute():
                # DePSI resolves relative names against its working directory;
                # the staged copy pins them to the original input-file directory
                # so the baseline run cannot silently read the wrong tree.
                path_i = (input_path.parent / path_i).resolve()
                relativized.append({'field': i, 'original': cells[i], 'staged': str(path_i)})
            if not path_i.is_file():
                raise ValueError(f'input list references a missing raster: {path_i}')
            cells[i] = str(path_i)
        input_lines.append(cells)
    out = Path(args.run_dir).resolve()
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'work'
    work.mkdir()
    (out / 'params.original.txt').write_bytes(params.read_bytes())
    (out / 'input.original.txt').write_bytes(input_path.read_bytes())
    for i, cells in enumerate(input_lines):
        if cells[3].lower() != 'dummy':
            cells[5] = str(work / f'scene_{i:04d}')
    staged_input = out / 'input.staged.txt'
    staged_input.write_text('\n'.join(' '.join(c) for c in input_lines) + '\n')
    staged_param = out / 'params.staged.txt'
    staged_text = re.sub(r"(?m)^\s*input_file\s*=.*$", "input_file = '" + str(staged_input).replace("'", "''") + "'", params.read_text())
    staged_param.write_text(staged_text)
    quote = lambda s: str(s).replace("'", "''")
    command = [args.matlab, '-batch', f"addpath('{quote(ROOT / 'scripts' / 'matlab')}'); capture_baseline('{quote(out / 'config.json')}')"]
    config = dict(schema=1, entry_file=str(entry), param_file=str(staged_param), work_dir=str(work),
                  result_file=str(out/'matlab-result.json'), profile_file=str(out/'profile.mat'),
                  profile=args.profile, entry_sha256=digest(entry), param_sha256=digest(params),
                  staged_param_sha256=digest(staged_param), input_sha256=digest(input_path),
                  resolved_relative_paths=relativized)
    repo = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    source_root = entry.parent.parent if entry.parent.name == 'main' else entry.parent
    source_manifest = {str(p.relative_to(source_root)): digest(p)
                       for p in sorted(source_root.rglob('*.m')) if '.git' not in p.parts and '.svn' not in p.parts}
    config['repository_revision'] = (getattr(repo, 'stdout', '') or '').strip() if repo.returncode == 0 else None
    config['source_manifest'] = source_manifest
    config['thread_environment'] = {k: os.environ.get(k) for k in
        ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS')}
    (out/'config.json').write_text(json.dumps(config, ensure_ascii=False, indent=2))
    report = dict(schema=1, status='prepared', executed=False, command=command, config=config, preflight=check)
    if args.execute:
        started = time.perf_counter()
        with (out/'matlab.log').open('w') as log:
            try:
                proc = subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
                code = proc.returncode
            except (subprocess.TimeoutExpired, OSError) as exc:
                code = -1
                report['error'] = str(exc)
        report.update(executed=True, returncode=code, elapsed_s=time.perf_counter()-started, status='failed')
        result = out/'matlab-result.json'
        if code == 0 and result.is_file():
            captured = json.loads(result.read_text())
            report['matlab_result'] = captured
            if captured.get('status') == 'completed':
                report['status'] = 'completed'
    (out/'run.json').write_text(json.dumps(report, ensure_ascii=False, indent=2))
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)

    def with_output(name, **kwargs):
        # --output must be accepted after the subcommand (argparse scoping)
        return sub.add_parser(name, **kwargs)

    env = with_output('environment')
    env.add_argument('--output', help='New JSON report path (never overwritten)')
    pre = with_output('preflight')
    pre.add_argument('params'); pre.add_argument('--map', help='JSON prefix mapping')
    pre.add_argument('--sample-pixels', type=int, default=0)
    pre.add_argument('--output', help='New JSON report path (never overwritten)')
    conv = with_output('from-mat', help='Convert explicit numeric fields from MATLAB -v7 MAT to NPZ')
    conv.add_argument('source'); conv.add_argument('target')
    conv.add_argument('--field', action='append', required=True)
    conv.add_argument('--output', help='New JSON report path (never overwritten)')
    comp = with_output('compare')
    comp.add_argument('left'); comp.add_argument('right')
    comp.add_argument('--output', help='New JSON report path (never overwritten)')
    comp.add_argument('--meta-left', required=True); comp.add_argument('--meta-right', required=True)
    comp.add_argument('--field', action='append', required=True)
    comp.add_argument('--atol', type=float, default=1e-8); comp.add_argument('--rtol', type=float, default=1e-7)
    mat = with_output('matlab')
    mat.add_argument('--output', help='New JSON report path (never overwritten)')
    mat.add_argument('--entry', required=True); mat.add_argument('--params', required=True)
    mat.add_argument('--run-dir', required=True); mat.add_argument('--matlab', default='matlab')
    mat.add_argument('--profile', action='store_true'); mat.add_argument('--execute', action='store_true')
    mat.add_argument('--timeout', type=float, default=3600)
    args = parser.parse_args(argv)
    try:
        if args.output and Path(args.output).exists():
            raise ValueError('Report exists; choose a new --output path before running an action.')
        if args.action == 'environment':
            report = environment()
        elif args.action == 'preflight':
            if not 0 <= args.sample_pixels <= 4096:
                raise ValueError('sample-pixels must be 0..4096 (bounded, not a full scan)')
            mapping = json.loads(Path(args.map).read_text()) if args.map else {}
            if not isinstance(mapping, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in mapping.items()):
                raise ValueError('path map must contain string prefixes')
            report = preflight(args.params, mapping, args.sample_pixels)
        elif args.action == 'from-mat':
            from scipy.io import loadmat
            target = Path(args.target)
            if target.exists():
                raise ValueError('Conversion target exists')
            if target.suffix != '.npz':
                raise ValueError('Conversion target must end in .npz')
            data = loadmat(args.source)
            names = ['ids'] + (['epochs'] if 'epochs' in data else []) + args.field
            arrays = {}
            for name in names:
                value = np.asarray(data[name])
                if value.dtype.kind not in 'biufc':
                    raise ValueError(f'{name}: numeric arrays only; use numeric acquisition dates YYYYMMDD')
                arrays[name] = value.ravel() if name == 'epochs' else value
            identities(arrays['ids'])
            with target.open('xb') as f:
                np.savez(f, **arrays)
            report = dict(schema=1, status='pass', scope='format conversion only',
                          source_sha256=digest(args.source), fields=names, target=str(target))
        elif args.action == 'compare':
            if args.atol < 0 or args.rtol < 0 or not np.isfinite([args.atol, args.rtol]).all():
                raise ValueError('tolerances must be finite and nonnegative')
            report = compare(args.left, args.right, args.meta_left, args.meta_right, args.field, args.atol, args.rtol)
        else:
            report = matlab_plan(args)
        emit(report, args.output)
        return 0 if report['status'] in ('pass', 'prepared', 'completed') else (1 if report['status'] in ('mismatch', 'failed') else 2)
    except (ValueError, OSError, KeyError, TypeError, IndexError) as exc:
        emit(dict(schema=1, status='blocked', error=str(exc)))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
