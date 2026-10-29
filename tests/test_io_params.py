"""Tests for the parameter-file parser (behaviour of readinput.m)."""

from __future__ import annotations

import numpy as np
import pytest

from depsipy.io.params import ParamError, Params, load_params, parse_param_file, parse_value


def _write(tmp_path, text, name="param_file.txt"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_scalars_ranges_and_arrays(tmp_path):
    p = _write(
        tmp_path,
        """
        % comment
        max_mem_buffer = 50e7
        processing_groups = [1:9]
        stc_min_max = [30,100]
        std_param = [30,0.005,1,0.02,0.01,0.01,0.005];
        empty_thing = []
        flag = 'yes'
        """,
    )
    d = parse_param_file(p)
    assert d["max_mem_buffer"] == pytest.approx(50e7)
    np.testing.assert_array_equal(d["processing_groups"], np.arange(1, 10))
    np.testing.assert_array_equal(d["stc_min_max"], [30, 100])
    assert d["std_param"].shape == (7,)
    assert d["std_param"][1] == pytest.approx(0.005)
    assert d["empty_thing"].size == 0
    assert d["flag"] == "yes"


def test_value_spanning_several_lines(tmp_path):
    p = _write(
        tmp_path,
        """
        std_param = [30,
                     0.005,
                     1]
        next_par = 3
        """,
    )
    d = parse_param_file(p)
    np.testing.assert_array_equal(d["std_param"], [30, 0.005, 1])
    assert d["next_par"] == 3


def test_range_with_step_and_expression(tmp_path):
    p = _write(tmp_path, "a = [1:2:9]\nb = 2*pi\nc = -1\n")
    d = parse_param_file(p)
    np.testing.assert_array_equal(d["a"], [1, 3, 5, 7, 9])
    assert d["b"] == pytest.approx(2 * np.pi)
    assert d["c"] == -1


def test_quoted_path_with_spaces_and_percent(tmp_path):
    p = _write(tmp_path, "input_file = '/data/50%_crop/input.txt'\n")
    d = parse_param_file(p)
    assert d["input_file"] == "/data/50%_crop/input.txt"


def test_duplicate_parameter_is_an_error(tmp_path):
    p = _write(tmp_path, "a = 1\na = 2\n")
    with pytest.raises(ParamError):
        parse_param_file(p)


def test_parse_value_direct():
    assert parse_value("42") == 42
    assert parse_value("42.5") == 42.5
    assert parse_value("'asc'") == "asc"
    assert parse_value("[]").size == 0
    np.testing.assert_array_equal(parse_value("[1,2,3]"), [1, 2, 3])


def test_params_typed_access(tmp_path):
    p = _write(tmp_path, "n = 3\ns = 'boot'\nv = [0.5,1.5]\n")
    params = load_params(p)
    assert params.int("n") == 3
    assert params.str("s") == "boot"
    npt = params.floats("v")
    np.testing.assert_allclose(npt, [0.5, 1.5])
    assert params.is_empty("nope")
    assert params.get("missing", 7) == 7
    with pytest.raises(ParamError):
        params.get("missing")


def test_real_depsi_param_file():
    """The parameter template shipped with the toolbox must parse."""
    from pathlib import Path

    ref = Path(__file__).resolve().parents[1] / "reference" / "templates" / "param_file.txt"
    if not ref.exists():  # pragma: no cover
        pytest.skip("reference template not present")
    d = parse_param_file(ref)
    assert d["ps_method"] == "boot"
    assert d["psc_threshold"] == pytest.approx(0.3)
    np.testing.assert_array_equal(d["std_param"], [30, 0.005, 1, 0.02, 0.01, 0.01, 0.005])
    np.testing.assert_array_equal(d["stc_min_max"], [30, 100])


def test_user_project_param_file():
    """The parameter file of the user's Spacety project must parse."""
    from pathlib import Path

    ref = (
        Path(__file__).resolve().parents[1]
        / "workspace"
        / "Doris干涉加时序depsi"
        / "depsi时序"
        / "param_file.txt"
    )
    if not ref.exists():  # pragma: no cover
        pytest.skip("user project parameter file not present")
    d = parse_param_file(ref)
    assert d["sensor"] == "spacety"
    assert d["ps_method"] == "boot"
    assert d["Ncon"] == 8
    np.testing.assert_array_equal(d["processing_groups"], np.arange(1, 10))
