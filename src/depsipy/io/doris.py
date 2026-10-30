"""Reader for Doris / RapidSAR products.

The interferometric chain in ``workspace/Doris干涉加时序depsi/`` (Doris + RapidSAR)
is what produces DePSI's inputs, and it writes them in its own conventions:

* a **MFF text header** (``IMAGE_FILE_FORMAT = MFF``, ``IMAGE_LINES``,
  ``LINE_SAMPLES``, ``FILE_TYPE``, ``BYTE_ORDER``) next to the data file;
* pixel-interleaved raw data (``.x00`` complex_real4, ``.r00`` real4,
  ``.j00`` complex_short), which is byte-identical to the DePSI raw layout;
* a **``.res`` result file** with ``key: value`` lines holding the radar
  parameters, the crop window and the precise orbit (``NUMBER_OF_DATAPOINTS``
  followed by ``[time x y z]`` rows in ECEF).

Parsing those three formats here is what makes the Doris products usable
directly (``depsipy.io.raw`` already reads the pixel data); the geocoding itself
lives in :mod:`depsipy.steps.geocoding`, mirroring ``metadataNew.m`` /
``intrp_orbit.m`` / ``lph2xyz.m``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .raw import RawFile, file_size, parse_format

__all__ = [
    "DorisMetadata",
    "DorisRaster",
    "MFFHeader",
    "doris_format",
    "read_doris_metadata",
    "read_mff_header",
    "read_res",
    "res_value",
]

#: Doris/RapidSAR data-file suffix -> DePSI format string (``SAR_write`` in
#: ``sar_image_crop.py`` uses the same three output formats)
_SUFFIX_FORMATS = {
    ".x00": "cpxfloat32",  # complex_real4
    ".r00": "float32",  # real4
    ".j00": "cpxint16",  # complex_short
}


def doris_format(path) -> str:
    """Return the DePSI format of a Doris data file from its suffix."""
    suffix = Path(path).suffix.lower()
    if suffix not in _SUFFIX_FORMATS:
        raise ValueError(
            f"unknown Doris data-file suffix {suffix!r}; expected one of "
            f"{', '.join(sorted(_SUFFIX_FORMATS))}"
        )
    return _SUFFIX_FORMATS[suffix]


# ---------------------------------------------------------------------------
# MFF header
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MFFHeader:
    """The subset of a Doris MFF header that matters for reading the raster."""

    path: Path
    lines: int
    samples: int
    file_type: str = "IMAGE"
    byte_order: str = "LSB"


def _header_path(path) -> Path:
    """Return the ``.hdr`` path next to ``path`` (or ``path`` itself)."""
    path = Path(path)
    if path.suffix.lower() == ".hdr":
        return path
    direct = path.with_suffix(".hdr")
    if direct.exists():
        return direct
    # Doris appends the header to the *data* name, not to the suffix:
    # ``mergedIfgBursts.raw`` -> ``mergedIfgBursts.hdr`` is the usual case.
    return direct


def read_mff_header(path) -> MFFHeader:
    """Parse a Doris MFF header (``IMAGE_FILE_FORMAT = MFF``)."""
    hdr = _header_path(path)
    if not hdr.exists():
        raise FileNotFoundError(f"MFF header not found for {path}: {hdr}")
    values: dict[str, str] = {}
    for raw in hdr.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().upper()] = value.strip()
    if "IMAGE_LINES" not in values or "LINE_SAMPLES" not in values:
        raise ValueError(f"{hdr}: not an MFF image header (IMAGE_LINES/LINE_SAMPLES missing)")
    return MFFHeader(
        path=hdr,
        lines=int(values["IMAGE_LINES"]),
        samples=int(values["LINE_SAMPLES"]),
        file_type=values.get("FILE_TYPE", "IMAGE"),
        byte_order=values.get("BYTE_ORDER", "LSB"),
    )


# ---------------------------------------------------------------------------
# Doris result file (``.res``)
# ---------------------------------------------------------------------------


_RES_KEY = re.compile(r"^\s*([^:]+?)\s*:\s*(.*)$", re.UNICODE)


def read_res(path) -> dict[str, str]:
    """Parse a Doris ``.res`` file into ``{key: value}`` (values kept as text).

    The MATLAB side tokenises on ``:`` and matches keys by *prefix*
    (``strmatch``); the same semantics are offered by :func:`res_value`.
    """
    path = Path(path)
    out: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _RES_KEY.match(raw)
        if not m:
            continue
        key, value = m.group(1).strip(), m.group(2).strip()
        if not key or key.startswith("*"):
            continue
        out.setdefault(key, value)
    return out


def res_value(res: dict[str, str], *prefixes: str, default=None):
    """Value of the first key starting with any of ``prefixes`` (``strmatch``)."""
    for prefix in prefixes:
        for key, value in res.items():
            if key.startswith(prefix):
                return value
    if default is not None:
        return default
    raise KeyError(f"none of {prefixes!r} found in the Doris result file")


# ---------------------------------------------------------------------------
# metadata container (port of ``metadataNew.m``)
# ---------------------------------------------------------------------------


def _res_float(res, *prefixes, default=None) -> float:
    value = res_value(res, *prefixes, default=None if default is None else str(default))
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"cannot parse {prefixes[0]!r} value {value!r} as a number") from exc


def _first_azimuth_time_seconds(text: str) -> float:
    """``19-OCT-1997 09:32:27.124`` -> seconds of day (``metadataNew.m``)."""
    m = re.search(r"(\d{1,2}):(\d{2}):(\d{2}(?:\.\d+)?)", text)
    if not m:
        raise ValueError(f"cannot parse the first azimuth time {text!r}")
    hh, mm, ss = (float(v) for v in m.groups())
    return hh * 3600.0 + mm * 60.0 + ss


@dataclass(frozen=True)
class DorisMetadata:
    """Radar parameters, crop window and orbit of one Doris result file.

    Field names follow the ``container`` vector of ``metadataNew.m``; the
    indices in the docstring are 1-based to make the mapping obvious.
    """

    path: Path
    orbit: np.ndarray  # (N, 4) [time, x, y, z] ECEF, seconds
    ta1: float  # container(1): seconds of day of the first azimuth line
    tr1: float  # container(2): one-way range time to the first pixel [s]
    prf: float  # container(3): PRF [Hz], multilook-compensated
    rsr: float  # container(4): one-way range sampling rate [Hz]
    wavelength: float  # container(5): radar wavelength [m]
    center_lat: float  # container(6): scene centre latitude [deg]
    center_lon: float  # container(7): scene centre longitude [deg]
    ovs_az: float  # container(8): 1/multilook factor, azimuth
    ovs_rg: float  # container(9): 1/multilook factor, range
    n_lines: float  # container(10): original lines, multilook-compensated
    n_pixels: float  # container(11): original pixels, multilook-compensated
    first_line: float  # container(12)
    last_line: float  # container(13)
    mid_line: float  # container(16): start point of the orbit interpolation
    first_pixel: float  # container(14)
    last_pixel: float  # container(15)
    dc_constant: float  # container(17) [Hz]
    dc_linear: float  # container(18) [Hz/s]
    dc_quadratic: float  # container(19) [Hz/s/s]
    data_file: str = ""  # containerChar(1)
    data_format: str = ""  # containerChar(2)
    raw: dict[str, str] = field(default_factory=dict, repr=False)


def read_doris_metadata(path) -> DorisMetadata:
    """Read the metadata container from a Doris ``.res`` file."""
    res = read_res(path)
    ovs_az = 1.0 / _res_float(res, "Multilookfactor_azimuth_direction", default=1.0)
    ovs_rg = 1.0 / _res_float(res, "Multilookfactor_range_direction", default=1.0)

    ta1 = _first_azimuth_time_seconds(res_value(res, "First_pixel_azimuth_time (UTC)"))
    tr1 = _res_float(res, "Range_time_to_first_pixel (2way) (ms)") / 2.0 / 1000.0
    prf = _res_float(
        res,
        "Pulse_Repetition_Frequency (computed, Hz)",
        "Pulse_Repetition_Frequency (actual, Hz)",
    )
    rsr = _res_float(
        res, "Range_sampling_rate (computed, MHz)", "Range_sampling_rate (leaderfile, MHz)"
    )
    wavelength = _res_float(res, "Radar_wavelength")
    center_lat = _res_float(res, "Scene_centre_latitude")
    center_lon = _res_float(res, "Scene_centre_longitude")

    first_line = _res_float(
        res, "First_line (w.r.t. original_image)", "First_line (w.r.t. original_master)"
    )
    last_line = _res_float(
        res, "Last_line (w.r.t. original_image)", "Last_line (w.r.t. original_master)"
    )
    first_pixel = _res_float(
        res, "First_pixel (w.r.t. original_image)", "First_pixel (w.r.t. original_master)"
    )
    last_pixel = _res_float(
        res, "Last_pixel (w.r.t. original_image)", "Last_pixel (w.r.t. original_master)"
    )
    n_lines = _res_float(res, "Number_of_lines_original")
    n_pixels = _res_float(res, "Number_of_pixels_original")

    # metadataNew.m converts the crop window to zero-based borders and accounts
    # for the multilook/oversampling factors.
    first_line = first_line * ovs_az - 1.0 if first_line > 1 else first_line
    last_line = last_line * ovs_az
    first_pixel = first_pixel * ovs_rg - 1.0 if first_pixel > 1 else first_pixel
    last_pixel = last_pixel * ovs_rg

    n_points = int(_res_float(res, "NUMBER_OF_DATAPOINTS"))
    orbit = _orbit_datapoints(path, n_points)

    return DorisMetadata(
        path=Path(path),
        orbit=orbit,
        ta1=ta1,
        tr1=tr1,
        prf=prf * ovs_az,
        rsr=rsr * 2.0 * ovs_rg * 1e6,
        wavelength=wavelength,
        center_lat=center_lat,
        center_lon=center_lon,
        ovs_az=ovs_az,
        ovs_rg=ovs_rg,
        n_lines=n_lines * ovs_az,
        n_pixels=n_pixels * ovs_rg,
        first_line=first_line,
        last_line=last_line,
        mid_line=(first_line + last_line) / 2.0,
        first_pixel=first_pixel,
        last_pixel=last_pixel,
        dc_constant=_res_float(res, "Xtrack_f_DC_constant (Hz, early edge)", default=0.0),
        dc_linear=_res_float(res, "Xtrack_f_DC_linear (Hz/s, early edge)", default=0.0),
        dc_quadratic=_res_float(res, "Xtrack_f_DC_quadratic (Hz/s/s, early edge)", default=0.0),
        data_file=str(res_value(res, "Data_output_file", default="")).strip(),
        data_format=str(res_value(res, "Data_output_format", default="")).strip(),
        raw=res,
    )


def _orbit_datapoints(path, n_points: int) -> np.ndarray:
    """Rows following ``NUMBER_OF_DATAPOINTS``: ``[time x y z]`` in ECEF."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for i, raw in enumerate(lines):
        if raw.split(":", 1)[0].strip().startswith("NUMBER_OF_DATAPOINTS"):
            start = i + 1
            break
    if start is None:
        raise ValueError(f"{path}: NUMBER_OF_DATAPOINTS not found")
    rows = []
    for raw in lines[start : start + n_points]:
        parts = raw.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: malformed orbit row {raw!r}")
        rows.append([float(v) for v in parts[:4]])
    if len(rows) != n_points:
        raise ValueError(f"{path}: expected {n_points} orbit rows, got {len(rows)}")
    orbit = np.asarray(rows, dtype=float)
    if np.any(np.diff(orbit[:, 0]) <= 0):
        raise ValueError(f"{path}: orbit epochs must be strictly increasing")
    return orbit


# ---------------------------------------------------------------------------
# raster
# ---------------------------------------------------------------------------


@dataclass
class DorisRaster:
    """A Doris raster described by its MFF header (1-based crop reads)."""

    path: Path
    lines: int
    samples: int
    fmt: str
    header: MFFHeader

    @classmethod
    def open(cls, path, *, fmt: str | None = None, header: MFFHeader | None = None) -> "DorisRaster":
        header = read_mff_header(path) if header is None else header
        fmt = doris_format(path) if fmt is None else fmt
        return cls(path=Path(path), lines=header.lines, samples=header.samples, fmt=fmt, header=header)

    def read(self, r0: int = 1, rN: int | None = None, c0: int = 1, cN: int | None = None) -> np.ndarray:
        """Read rows ``r0..rN`` / columns ``c0..cN`` (1-based, inclusive)."""
        rf = parse_format(self.fmt)
        expected = self.lines * self.samples * rf.bytes_per_pixel
        size = file_size(self.path)
        if size != expected:
            raise ValueError(
                f"{self.path}: {size} bytes on disk, MFF header implies {expected} "
                f"({self.lines} lines x {self.samples} samples, {self.fmt})"
            )
        width = int(round(size / (self.lines * rf.bytes_per_pixel)))
        with RawFile(self.path, self.lines, self.fmt) as f:
            if f.width != self.samples or width != self.samples:
                raise ValueError(
                    f"{self.path}: inferred width {f.width} does not match "
                    f"MFF LINE_SAMPLES {self.samples}"
                )
            return f.read(r0, rN, c0, cN)
