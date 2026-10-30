"""Doris product ingestion and radar geocoding tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from depsipy.io.doris import (
    DorisRaster,
    doris_format,
    read_doris_metadata,
    read_mff_header,
    read_res,
)
from depsipy.io.raw import write_raw
from depsipy.steps.geocoding import (
    latlonh_to_xyz,
    lph2xyz,
    orbit_state,
    radar_to_geodetic,
    xyz2lp,
)

A = 6378137.0
REFERENCE_RES = Path(__file__).resolve().parents[1] / "reference/matlab/geocoding_v0.8/example.res"


def _write_res(path: Path, *, n_orbit: int = 5) -> Path:
    lines = [
        " MASTER RESULTFILE:  master.res",
        "Start_process_control",
        "crop:                   1",
        "End_process_control",
        "*_Start_readfiles:",
        "Radar_wavelength (m):                           0.05656461471698114",
        "First_pixel_azimuth_time (UTC):                 01-JAN-2020 00:00:00.000",
        "Pulse_Repetition_Frequency (computed, Hz):      3359.806849229288",
        "Range_time_to_first_pixel (2way) (ms):          5.55976",
        "Range_sampling_rate (computed, MHz):            37.92491979726448",
        "Scene_centre_latitude:                          49.5750389",
        "Scene_centre_longitude:                         20.2600899",
        "Xtrack_f_DC_constant (Hz, early edge):          197.9807281",
        "Xtrack_f_DC_linear (Hz/s, early edge):          62338.3984375",
        "Xtrack_f_DC_quadratic (Hz/s/s, early edge):     -485652224.00000",
        "Number_of_lines_original:                       26441",
        "Number_of_pixels_original:                      4904",
        "Data_output_file:                               master_crop.raw",
        "Data_output_format:                             complex_short",
        "First_line (w.r.t. original_image):             131",
        "Last_line (w.r.t. original_image):              10130",
        "First_pixel (w.r.t. original_image):            6023",
        "Last_pixel (w.r.t. original_image):             8022",
        f"NUMBER_OF_DATAPOINTS:                           {n_orbit}",
    ]
    for k in range(n_orbit):
        lines.append(f"{34337.0 + k:.6f}    {4198191.8 + k:.3f}    {1931211.741 + k:.3f}    5465346.293")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def test_read_res_and_metadata_synthetic(tmp_path):
    res_path = _write_res(tmp_path / "master.res")
    res = read_res(res_path)
    assert res["Radar_wavelength (m)"] == "0.05656461471698114"
    assert res["Data_output_format"] == "complex_short"

    meta = read_doris_metadata(res_path)
    assert meta.ta1 == pytest.approx(0.0)  # 00:00:00.000
    assert meta.prf == pytest.approx(3359.806849229288)
    assert meta.tr1 == pytest.approx(5.55976 / 2 / 1000)
    assert meta.rsr == pytest.approx(37.92491979726448 * 2 * 1e6)
    assert meta.wavelength == pytest.approx(0.05656461471698114)
    assert (meta.center_lat, meta.center_lon) == pytest.approx((49.5750389, 20.2600899))
    # metadataNew.m keeps the zero-based crop border (first_line - 1)
    assert meta.first_line == pytest.approx(130.0)
    assert meta.mid_line == pytest.approx((130.0 + 10130.0) / 2.0)
    assert meta.orbit.shape == (5, 4)
    assert np.all(np.diff(meta.orbit[:, 0]) > 0)
    assert meta.data_file == "master_crop.raw"
    assert meta.data_format == "complex_short"


@pytest.mark.skipif(not REFERENCE_RES.exists(), reason="reference/ is local-only")
def test_read_res_reference_example():
    meta = read_doris_metadata(REFERENCE_RES)
    assert meta.ta1 == pytest.approx(34347.124)  # 09:32:27.124
    assert meta.prf == pytest.approx(3359.806849229288)
    assert meta.wavelength == pytest.approx(0.05656461471698114)
    assert meta.orbit.shape == (37, 4)
    assert meta.first_line == pytest.approx(130.0)
    assert meta.last_line == pytest.approx(10130.0)
    assert meta.mid_line == pytest.approx(5130.0)


def test_mff_header_and_raster(tmp_path):
    data = np.arange(30, dtype=np.float32).reshape(6, 5)
    raw = tmp_path / "ifg.r00"
    write_raw(raw, data, "float32")
    (tmp_path / "ifg.hdr").write_text(
        "IMAGE_FILE_FORMAT = MFF\nFILE_TYPE = IMAGE\nIMAGE_LINES = 6\n"
        "LINE_SAMPLES = 5\nBYTE_ORDER = LSB\nEND\n",
        encoding="utf-8",
    )
    header = read_mff_header(raw)
    assert (header.lines, header.samples) == (6, 5)

    raster = DorisRaster.open(raw)
    block = raster.read(2, 4, 2, 3)
    np.testing.assert_array_equal(block, data[1:4, 1:3])
    assert doris_format(raw) == "float32"
    assert doris_format(tmp_path / "ifg.x00") == "cpxfloat32"


def test_mff_size_mismatch_is_rejected(tmp_path):
    data = np.zeros((5, 5), dtype=np.float32)
    raw = tmp_path / "bad.r00"
    write_raw(raw, data, "float32")
    (tmp_path / "bad.hdr").write_text(
        "IMAGE_FILE_FORMAT = MFF\nIMAGE_LINES = 6\nLINE_SAMPLES = 5\nEND\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="bytes on disk"):
        DorisRaster.open(raw).read()


def _straight_orbit_meta(center_lat: float = 0.05):
    from depsipy.io.doris import DorisMetadata

    v = 7000.0  # m/s
    t = np.arange(0.0, 31.0, 1.0)
    orbit = np.column_stack([t, np.full_like(t, A + 700e3), v * t, np.zeros_like(t)])
    return DorisMetadata(
        path=Path("<synthetic>"),
        orbit=orbit,
        ta1=0.0,
        tr1=2.0e-3,
        prf=3000.0,
        rsr=30e6,
        wavelength=0.0555,
        center_lat=center_lat,
        center_lon=0.0,
        ovs_az=1.0,
        ovs_rg=1.0,
        n_lines=10000.0,
        n_pixels=20000.0,
        first_line=1.0,
        last_line=10000.0,
        mid_line=5000.0,
        first_pixel=1.0,
        last_pixel=20000.0,
        dc_constant=0.0,
        dc_linear=0.0,
        dc_quadratic=0.0,
    )


def test_geocoding_roundtrip():
    meta = _straight_orbit_meta()
    lat, lon, height = 0.20, 0.01, 0.0
    target = latlonh_to_xyz(lat, lon, height)

    line, pixel = xyz2lp(target, meta)
    _, orbit_line = orbit_state(meta, line)
    xyz = lph2xyz(line, pixel, [height], meta, orbit_line)
    assert np.linalg.norm(xyz[0] - target) < 1e-4

    geo = radar_to_geodetic(line, pixel, height, meta)[0]
    assert geo[0] == pytest.approx(np.rad2deg(lat), abs=1e-6)
    assert geo[1] == pytest.approx(np.rad2deg(lon), abs=1e-6)
    assert geo[2] == pytest.approx(height, abs=1e-3)


def test_geocoding_roundtrip_with_height():
    meta = _straight_orbit_meta()
    lat, lon, height = 0.20, 0.01, 500.0
    target = latlonh_to_xyz(lat, lon, height)
    line, pixel = xyz2lp(target, meta)
    _, orbit_line = orbit_state(meta, line)
    xyz = lph2xyz(line, pixel, [height], meta, orbit_line)
    assert np.linalg.norm(xyz[0] - target) < 1e-3


def test_orbit_interpolation_matches_a_straight_line():
    meta = _straight_orbit_meta()
    _, state = orbit_state(meta, [1000.0, 5000.0, 9000.0])
    # parabolic fit of a straight line is exact: constant velocity, zero acceleration
    np.testing.assert_allclose(state[:, 4], 0.0, atol=1e-6)
    np.testing.assert_allclose(state[:, 5], 7000.0, atol=1e-6)
    np.testing.assert_allclose(state[:, 7:10], 0.0, atol=1e-9)
    np.testing.assert_allclose(state[:, 0], (np.array([1000.0, 5000.0, 9000.0]) - 1) / 3000.0)
