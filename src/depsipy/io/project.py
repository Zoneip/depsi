"""Project description: the ``input_file`` and its derived geometry.

Port of the "Read input" block of ``depsi.m`` (columns, master handling, sorting,
`crop_final`) plus the sensor constants.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["SENSORS", "InputStack", "read_input_file", "sensor_constants"]

#: (wavelength [m], satellite velocity [m/s]) for the sensors ``depsi.m`` knows
SENSORS = {
    "ers": (0.0565646, 7550.0),
    "asar": (0.0562356, 7550.0),
    "rsat1": (0.0564105, 7550.0),
    "rsat2": (0.055465772433, 7550.0),
    "tsx": (0.031000, 7600.0),
    "sentinel-1": (0.031000, 7590.0),
    "lt": (0.031000, 7590.0),
    "spacety": (0.055517120823757024, 7690.0),
}


def sensor_constants(sensor: str) -> tuple[float, float, float]:
    """Return ``(lambda, satellite_velocity, m2ph)`` for ``sensor``."""
    key = str(sensor).strip().lower()
    if key not in SENSORS:
        raise ValueError(
            f"unsupported sensor {sensor!r}; known: {', '.join(sorted(SENSORS))}"
        )
    lam, vel = SENSORS[key]
    return lam, vel, -4 * np.pi / lam


@dataclass
class InputStack:
    """Contents of an ``input_file`` in DePSI order (master removed)."""

    orbit: np.ndarray  # (Nifgs,) str
    date: np.ndarray  # (Nifgs,) str
    slc: list[str]
    ifg: list[str]
    h2ph: list[str]
    out: list[str]
    Bperp: np.ndarray
    Btemp: np.ndarray  # years
    Bdop: np.ndarray
    crop_in: np.ndarray  # (Nslc, 4) 1-based [l0 lN p0 pN]; row Nifgs = master
    crop_final: tuple[int, int, int, int]
    master_index: int  # 1-based
    slc_master: str = ""
    root: Path = field(default_factory=Path)

    def slc_paths(self) -> list[str]:
        """SLC filenames in DePSI order (``Nifgs`` interferometric scenes + master)."""
        return [self.resolve(f) for f in self.slc]

    @property
    def Nifgs(self) -> int:
        return self.Btemp.size

    @property
    def Nslc(self) -> int:
        return self.Nifgs + 1

    @property
    def Nslc_files(self) -> int:
        return len(self.slc)

    def resolve(self, name: str) -> str:
        """Absolute path of a file listed in the input file (relative to its dir)."""
        p = Path(name)
        return str(p if p.is_absolute() else self.root / p)


def read_input_file(path) -> InputStack:
    """Parse a DePSI ``input_file``.

    Columns: ``orbit date slc ifg h2ph out Bperp Btemp Bdop l0 lN p0 pN``.
    The master row carries ``dummy`` in the ``ifg``/``h2ph``/``out`` columns and is
    moved to the end, exactly as ``depsi.m`` does.
    """
    path = Path(path)
    rows: list[list[str]] = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) != 13:
            raise ValueError(
                f"{path}:{lineno}: expected 13 whitespace separated columns, got {len(parts)}"
            )
        rows.append(parts)
    if not rows:
        raise ValueError(f"{path}: empty input file")

    orbit = np.array([r[0] for r in rows])
    date = np.array([r[1] for r in rows])
    slc = [r[2] for r in rows]
    ifg = [r[3] for r in rows]
    h2ph = [r[4] for r in rows]
    out = [r[5] for r in rows]
    Bperp = np.array([float(r[6]) for r in rows])
    Btemp = np.array([float(r[7]) for r in rows])
    Bdop = np.array([float(r[8]) for r in rows])
    crop_in = np.array([[int(r[9]), int(r[10]), int(r[11]), int(r[12])] for r in rows], dtype=int)

    is_master = np.array([g.lower() == "dummy" for g in ifg])
    if is_master.sum() != 1:
        raise ValueError(
            f"{path}: exactly one row must be the master (ifg == 'dummy'), found {is_master.sum()}"
        )
    m = int(np.nonzero(is_master)[0][0])
    slc_master = slc[m]

    keep = np.arange(len(rows)) != m
    order = np.argsort(Btemp[keep], kind="stable")
    keep = keep.nonzero()[0][order]

    Btemp = Btemp[keep]
    if np.abs(Btemp).max() > 33:  # MATLAB: if Btemp is in days, convert to years
        Btemp = Btemp / 365.0

    # master index in the sorted order: last scene with Btemp < 0, +1
    neg = np.nonzero(Btemp < 0)[0]
    master_index = int(neg[-1]) + 1 if neg.size else 1

    crop_final = (
        int(crop_in[:, 0].max()),
        int(crop_in[:, 1].min()),
        int(crop_in[:, 2].max()),
        int(crop_in[:, 3].min()),
    )
    # MATLAB permutes only the first Nifgs rows of crop_in/filenames_slc and keeps
    # the master row in place, so crop_in has Nslc rows.
    crop_in_sorted = np.vstack([crop_in[keep], crop_in[m : m + 1]])
    return InputStack(
        orbit=orbit[keep],
        date=date[keep],
        slc=[slc[i] for i in keep] + [slc_master],
        ifg=[ifg[i] for i in keep],
        h2ph=[h2ph[i] for i in keep],
        out=[out[i] for i in keep],
        Bperp=Bperp[keep],
        Btemp=Btemp,
        Bdop=Bdop[keep],
        crop_in=crop_in_sorted,
        crop_final=crop_final,
        master_index=master_index,
        slc_master=slc_master,
        root=path.parent,
    )
