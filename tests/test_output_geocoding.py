"""Tests for output and geocoding basics (F4/W5)."""

from __future__ import annotations

import json

import numpy as np

from depsipy.io.raw import read_raw
from depsipy.steps.geocoding import latlonh_to_xyz, xyz_to_latlonh
from depsipy.steps.output import write_point_csv, write_point_stack, write_resfile, write_summary


def test_write_point_stack_roundtrip(tmp_path):
    az = np.array([1.0, 2.0, 3.0])
    r = np.array([4.0, 5.0, 6.0])
    values = np.array([[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]])
    path = write_point_stack(tmp_path / "points.raw", az, r, values)
    data = read_raw(path, 3, "float64")
    np.testing.assert_allclose(data, np.column_stack([az, r, values]))


def test_write_csv_and_resfile(tmp_path):
    az = np.array([1, 2])
    r = np.array([3, 4])
    values = np.array([[5.0], [6.0]])
    csv = write_point_csv(tmp_path / "points.csv", az, r, values, columns=["trend"])
    lines = csv.read_text().splitlines()
    assert lines[0] == "azimuth,range,trend"
    assert list(map(float, lines[1].split(","))) == [1.0, 3.0, 5.0]
    res = write_resfile(tmp_path / "resfile.txt", {"Npsc": 2, "Narcs": 1})
    assert "Npsc: 2" in res.read_text()
    summary = write_summary(tmp_path / "summary.json", {"status": "pass"})
    assert json.loads(summary.read_text())["status"] == "pass"


def test_geodetic_roundtrip():
    lat = np.array([0.0, 0.5, -1.0])
    lon = np.array([0.0, 1.2, -0.7])
    h = np.array([0.0, 1000.0, -30.0])
    xyz = latlonh_to_xyz(lat, lon, h)
    np.testing.assert_allclose(xyz[0], [6378137.0, 0.0, 0.0], atol=1e-6)
    back = xyz_to_latlonh(xyz)
    np.testing.assert_allclose(back[:, 0], lat, atol=1e-10)
    np.testing.assert_allclose(back[:, 1], lon, atol=1e-10)
    np.testing.assert_allclose(back[:, 2], h, atol=1e-6)
