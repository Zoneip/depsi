"""Tests for the raw binary I/O layer (byte-compatibility with freadbk/fwritebk)."""

from __future__ import annotations

import numpy as np
import pytest

from depsipy.io import raw


def test_roundtrip_float32(tmp_path):
    path = tmp_path / "a.raw"
    data = np.arange(12, dtype=np.float32).reshape(3, 4)
    n = raw.write_raw(path, data, "float32")
    assert n == 12
    back = raw.read_raw(path, 3, "float32")
    np.testing.assert_array_equal(back, data)


def test_roundtrip_complex_matches_matlab_layout(tmp_path):
    """fwritebk writes real/imag interleaved per pixel in row-major order."""
    path = tmp_path / "c.raw"
    data = (np.arange(6) + 1j * np.arange(6, 12)).astype(np.complex64).reshape(2, 3)
    raw.write_raw(path, data, "cpxfloat32")
    # explicit byte check: pixel (0,0) then (0,1) ... row by row
    flat = np.fromfile(path, dtype=np.float32)
    expected = np.empty(12, dtype=np.float32)
    expected[0::2] = data.real.ravel()
    expected[1::2] = data.imag.ravel()
    np.testing.assert_allclose(flat, expected)
    np.testing.assert_array_equal(raw.read_raw(path, 2, "cpxfloat32"), data)


def test_cropped_read_is_inclusive_and_one_based(tmp_path):
    path = tmp_path / "b.raw"
    data = np.arange(100, dtype=np.float64).reshape(10, 10)
    raw.write_raw(path, data, "float64")
    got = raw.read_raw(path, 10, "float64", r0=2, rN=4, c0=3, cN=5)
    np.testing.assert_array_equal(got, data[1:4, 2:5])


def test_default_crop_reads_everything(tmp_path):
    path = tmp_path / "c2.raw"
    data = np.random.randn(7, 5).astype(np.float32)
    raw.write_raw(path, data)
    np.testing.assert_array_equal(raw.read_raw(path, 7), data)


def test_infer_shape_detects_wrong_lines(tmp_path):
    path = tmp_path / "d.raw"
    raw.write_raw(path, np.zeros((4, 6), dtype=np.float32))
    lines, width = raw.infer_shape(path, 4, "float32")
    assert (lines, width) == (4, 6)
    with pytest.raises(ValueError):
        raw.infer_shape(path, 5, "float32")


def test_complex_int16_format(tmp_path):
    path = tmp_path / "e.raw"
    re = np.array([[1, -2], [3, -4]], dtype=np.int16)
    im = np.array([[5, -6], [7, -8]], dtype=np.int16)
    inter = np.empty((2, 4), dtype=np.int16)
    inter[:, 0::2] = re
    inter[:, 1::2] = im
    inter.tofile(path)
    got = raw.read_raw(path, 2, "cpxint16")
    np.testing.assert_array_equal(got, re + 1j * im)


def test_read_pixels_matches_full_read(tmp_path):
    path = tmp_path / "f.raw"
    data = np.random.randn(50, 40).astype(np.float32)
    raw.write_raw(path, data)
    az = np.array([1, 50, 17, 33])
    r = np.array([1, 40, 5, 22])
    got = raw.read_pixels(path, 50, "float32", az, r)
    np.testing.assert_array_equal(got, data[az - 1, r - 1])


def test_rawfile_mapping_and_close(tmp_path):
    path = tmp_path / "g.raw"
    data = np.random.randn(20, 30).astype(np.float32)
    raw.write_raw(path, data)
    with raw.RawFile(path, 20, "float32") as f:
        assert (f.lines, f.width) == (20, 30)
        np.testing.assert_array_equal(f.read(5, 8, 2, 4), data[4:8, 1:4])
        np.testing.assert_array_equal(f.row(3), data[2])
    assert f._mm is None


def test_iter_blocks_covers_all():
    assert list(raw.iter_blocks(10, 4)) == [(1, 4), (5, 8), (9, 10)]
    assert list(raw.iter_blocks(0, 4)) == []


def test_complex_width_counts_pixels_not_reals(tmp_path):
    """freadbk computes the width from bytes/elem where cpx counts 8 bytes."""
    path = tmp_path / "h.raw"
    raw.write_raw(path, np.ones((5, 7), dtype=np.complex64), "cpxfloat32")
    lines, width = raw.infer_shape(path, 5, "cpxfloat32")
    assert (lines, width) == (5, 7)
