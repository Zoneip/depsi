"""Result output helpers (step7): point tables, raw stacks, resfile, summary.

DePSI's original products are raw rasters plus ``_resfile.txt``; this module
keeps the same conventions for point time series and adds CSV/JSON for
interchange.  Higher-level raster/GeoTIFF output is intentionally left to a
later optional extra.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..io.raw import write_raw

__all__ = ["write_point_stack", "write_point_csv", "write_resfile", "write_summary"]


def write_point_stack(
    path,
    az: np.ndarray,
    r: np.ndarray,
    values: np.ndarray,
    *,
    fmt: str = "float64",
) -> Path:
    """Write ``[az r value_1 ... value_k]`` per point as a raw float stack."""
    az = np.asarray(az, dtype=float).ravel()
    r = np.asarray(r, dtype=float).ravel()
    val = np.atleast_2d(np.asarray(values, dtype=float))
    if val.shape[0] != az.size or az.size != r.size:
        raise ValueError("az, r and value rows must have the same size")
    data = np.column_stack([az, r, val])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_raw(path, data, fmt)
    return path


def write_point_csv(
    path,
    az: np.ndarray,
    r: np.ndarray,
    values: np.ndarray,
    *,
    columns: list[str] | None = None,
) -> Path:
    """Write a human-readable point table with ``azimuth,range`` as first columns."""
    az = np.asarray(az, dtype=float).ravel()
    r = np.asarray(r, dtype=float).ravel()
    val = np.atleast_2d(np.asarray(values, dtype=float))
    if val.shape[0] != az.size or az.size != r.size:
        raise ValueError("az, r and value rows must have the same size")
    names = list(columns or [f"v{i + 1}" for i in range(val.shape[1])])
    if len(names) != val.shape[1]:
        raise ValueError("number of column names must match values")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ",".join(["azimuth", "range", *names])
    np.savetxt(path, np.column_stack([az, r, val]), delimiter=",", header=header, comments="")
    return path


def write_resfile(path, entries: dict[str, object]) -> Path:
    """Write a DePSI-style ``_resfile.txt`` (key: value lines)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for key, value in entries.items():
            f.write(f"{key}: {value}\n")
    return path


def write_summary(path, summary: dict[str, object]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path
