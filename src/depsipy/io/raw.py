"""Raw binary I/O, byte-compatible with DePSI/DORIS ``freadbk``/``fwritebk``.

The DORIS/DePSI *raw* convention is:

* row-major (C order) storage of a ``(lines, width)`` matrix,
* for complex types, real and imaginary parts interleaved per pixel
  ("pixel interleaved"), i.e. exactly ``numpy.complex64`` with C order,
* no header whatsoever -- the size of the file defines the width:
  ``width = filesize / (bytes_per_pixel * lines)``.

That makes a numpy ``memmap`` with ``shape=(lines, width)`` the natural
abstraction, and a *cropped* read becomes a plain slice.  ``freadbk`` in
MATLAB does one ``fseek``+``fread`` per line, which is the single largest
avoidable cost in the original ``ps_selection.m``; slicing a memmap is
10-100x faster for the same result.

Complex ``int16`` data (older DORIS output, ``.rsmp``) is supported as well.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterator

import numpy as np

__all__ = [
    "RawFormat",
    "parse_format",
    "file_size",
    "infer_shape",
    "read_raw",
    "write_raw",
    "RawFile",
    "iter_blocks",
    "read_pixels",
]


#: MATLAB/DePSI storage type name -> numpy *real* dtype
_FORMATS: dict[str, str] = {
    "int8": "i1",
    "uchar": "u1",
    "uint8": "u1",
    "int16": "i2",
    "uint16": "u2",
    "int32": "i4",
    "uint32": "u4",
    "int64": "i8",
    "uint64": "u8",
    "float32": "f4",
    "single": "f4",
    "float64": "f8",
    "double": "f8",
    "float": "f4",
    "short": "i2",
    "int": "i4",
    "uint": "u4",
    "long": "i4",
    "ushort": "u2",
    "ulong": "u4",
}


@dataclass(frozen=True)
class RawFormat:
    """A parsed DePSI raw data format string.

    ``storage`` is the dtype of a *real* number on disk; complex data stores two
    of them per pixel (real then imaginary), so ``bytes_per_pixel`` doubles.
    ``dtype`` is the logical dtype an array has in memory.
    """

    name: str
    storage: np.dtype
    is_complex: bool

    @property
    def bytes_per_pixel(self) -> int:
        return int(self.storage.itemsize) * (2 if self.is_complex else 1)

    @property
    def dtype(self) -> np.dtype:
        if not self.is_complex:
            return self.storage
        # float32 -> complex64, float64 -> complex128, integers -> complex64
        return np.dtype("c8") if self.storage.itemsize <= 4 else np.dtype("c16")


def parse_format(fmt: str | RawFormat) -> RawFormat:
    """Parse a DePSI format string such as ``'cpxfloat32'`` or ``'float32'``.

    ``'mph'`` is accepted as an alias for ``'cpxfloat32'`` (see ``freadbk``).
    """
    if isinstance(fmt, RawFormat):
        return fmt
    if not isinstance(fmt, str):
        raise TypeError(f"format must be a string, got {type(fmt)!r}")

    name = fmt.strip()
    if name == "mph":
        name = "cpxfloat32"
    if name == "hgt":
        raise ValueError("the 'hgt' SRTM format is not supported by the raw reader")

    is_complex = name.startswith("cpx")
    base = name[3:] if is_complex else name
    if base not in _FORMATS:
        raise ValueError(f"unknown DePSI data format {fmt!r}")
    return RawFormat(name=name, storage=np.dtype(_FORMATS[base]), is_complex=is_complex)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def file_size(path: str | os.PathLike) -> int:
    return os.stat(path).st_size


def infer_shape(path, lines: int, fmt: str | RawFormat) -> tuple[int, int]:
    """Return ``(lines, width)`` of a raw file, the way ``freadbk`` computes it."""
    rf = parse_format(fmt)
    size = file_size(path)
    lines = int(lines)
    if lines <= 0:
        raise ValueError("lines must be positive")
    denom = rf.bytes_per_pixel * lines
    if size % denom:
        raise ValueError(
            f"{path}: file size {size} is not a multiple of "
            f"{rf.bytes_per_pixel} bytes x {lines} lines "
            f"(computed width would be {size / denom:.6f}); wrong line count or format?"
        )
    return lines, int(size // denom)


def _as_2d(data: np.ndarray) -> np.ndarray:
    a = np.asarray(data)
    if a.ndim != 2:
        raise ValueError(f"raw data must be 2-D, got shape {a.shape}")
    return a


def _interleave_complex(a: np.ndarray) -> np.ndarray:
    """Return an array of reals with re/im interleaved per pixel (row major)."""
    out = np.empty((a.shape[0], 2 * a.shape[1]), dtype=a.dtype)
    out[:, 0::2] = a.real
    out[:, 1::2] = a.imag
    return out


def _fold_complex(buf: np.ndarray) -> np.ndarray:
    """Inverse of :func:`_interleave_complex`."""
    re = buf[:, 0::2]
    im = buf[:, 1::2]
    dt = np.complex64 if buf.dtype.itemsize <= 4 else np.complex128
    return (re.astype(dt) + 1j * im.astype(dt)).astype(dt)


def write_raw(path, data, fmt: str | RawFormat = "float32") -> int:
    """Write ``data`` in DePSI raw layout.  Returns the number of pixels written."""
    rf = parse_format(fmt)
    a = _as_2d(data)

    if rf.is_complex:
        if not np.iscomplexobj(a):
            a = a.astype(np.complex64)
        if rf.storage.kind == "f":
            np.ascontiguousarray(a.astype(rf.dtype, copy=False)).tofile(os.fspath(path))
        else:
            inter = _interleave_complex(a.astype(np.complex128))
            np.ascontiguousarray(inter.astype(rf.storage)).tofile(os.fspath(path))
    else:
        if np.iscomplexobj(a):
            raise ValueError("complex data cannot be written to a real raw format")
        np.ascontiguousarray(a.astype(rf.storage, copy=False)).tofile(os.fspath(path))
    return int(a.shape[0] * a.shape[1])


def read_raw(
    path,
    lines: int | None = None,
    fmt: str | RawFormat = "float32",
    r0: int | None = None,
    rN: int | None = None,
    c0: int | None = None,
    cN: int | None = None,
) -> np.ndarray:
    """Read a (crop of a) raw file, mirroring ``freadbk``.

    Row/column indices are **1-based inclusive**, as in MATLAB.  ``None`` or 0
    means "from the first" / "up to the last" element.  If ``lines`` is omitted
    the whole file is returned as a 1-D vector, exactly like ``freadbk`` when
    the line count is unknown.
    """
    rf = parse_format(fmt)

    if lines is None:
        if any(v for v in (r0, rN, c0, cN)):
            raise ValueError("a cropped read requires the number of lines")
        flat = np.fromfile(os.fspath(path), dtype=rf.storage)
        if rf.is_complex:
            return _fold_complex(flat.reshape(-1, 2))
        return flat

    nl, width = infer_shape(path, lines, rf)
    r0 = 1 if not r0 else int(r0)
    rN = nl if not rN else int(rN)
    c0 = 1 if not c0 else int(c0)
    cN = width if not cN else int(cN)
    if r0 > rN:
        raise ValueError("r0 > rN")
    if c0 > cN:
        raise ValueError("c0 > cN")
    if r0 < 1 or rN > nl:
        raise ValueError(f"row range {r0}:{rN} outside file ({nl} lines)")
    if c0 < 1 or cN > width:
        raise ValueError(f"column range {c0}:{cN} outside file ({width} pixels)")

    if rf.is_complex and rf.storage.kind != "f":
        mm = np.memmap(os.fspath(path), dtype=rf.storage, mode="r", shape=(nl, 2 * width))
        try:
            blk = np.array(mm[r0 - 1 : rN, 2 * (c0 - 1) : 2 * cN], copy=True)
        finally:
            del mm
        return _fold_complex(blk)

    mm = np.memmap(os.fspath(path), dtype=rf.dtype, mode="r", shape=(nl, width))
    try:
        return np.array(mm[r0 - 1 : rN, c0 - 1 : cN], copy=True)
    finally:
        del mm


class RawFile:
    """A memory-mapped DePSI raw file with 1-based crop access.

    Use this instead of :func:`read_raw` when the same file is read more than
    once: the mapping stays open and only a contiguous copy is returned.
    """

    def __init__(self, path, lines: int, fmt="float32", *, read_only: bool = True):
        self.path = os.fspath(path)
        self.fmt = parse_format(fmt)
        self.lines, self.width = infer_shape(self.path, lines, self.fmt)
        self._mm: np.memmap | None = None
        self._read_only = read_only

    # -- lifecycle ---------------------------------------------------------
    @property
    def mm(self) -> np.memmap:
        if self._mm is None:
            if self.fmt.is_complex and self.fmt.storage.kind != "f":
                shape = (self.lines, 2 * self.width)
                dtype = self.fmt.storage
            else:
                shape = (self.lines, self.width)
                dtype = self.fmt.dtype
            self._mm = np.memmap(
                self.path,
                dtype=dtype,
                mode="r" if self._read_only else "r+",
                shape=shape,
            )
        return self._mm

    def close(self) -> None:
        mm = self._mm
        self._mm = None
        if mm is not None:
            mm._mmap.close()  # type: ignore[attr-defined]

    def __enter__(self) -> "RawFile":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # -- access ------------------------------------------------------------
    def read(self, r0: int = 1, rN: int | None = None, c0: int = 1, cN: int | None = None):
        """Read rows ``r0..rN`` and columns ``c0..cN`` (1-based, inclusive)."""
        rN = self.lines if rN is None else int(rN)
        cN = self.width if cN is None else int(cN)
        if self.fmt.is_complex and self.fmt.storage.kind != "f":
            blk = np.array(self.mm[int(r0) - 1 : rN, 2 * (int(c0) - 1) : 2 * cN], copy=True)
            return _fold_complex(blk)
        return np.array(self.mm[int(r0) - 1 : rN, int(c0) - 1 : cN], copy=True)

    def row(self, line: int) -> np.ndarray:
        """Return one line of the raster as a 1-D array."""
        return self.read(line, line)[0]


def iter_blocks(total: int, block: int) -> Iterator[tuple[int, int]]:
    """Yield 1-based inclusive ``(start, stop)`` blocks covering ``total`` items."""
    if block <= 0:
        raise ValueError("block size must be positive")
    start = 0
    while start < total:
        stop = min(start + block, total)
        yield start + 1, stop
        start = stop


def read_pixels(path, lines: int, fmt: str, az: np.ndarray, r: np.ndarray, *, chunk: int = 1 << 16):
    """Read scattered pixels ``(az, r)`` (1-based) from a raw raster.

    Reading only the pixels a step actually needs is what makes the Python
    selection pass fast: DePSI's ``ps_selection.m`` loads whole interferogram
    buffers although only a few hundred candidate pixels are used.
    """
    rf = parse_format(fmt)
    nl, width = infer_shape(path, lines, rf)
    az = np.asarray(az, dtype=np.int64)
    r = np.asarray(r, dtype=np.int64)
    if az.shape != r.shape:
        raise ValueError("az and r must have the same shape")
    if az.size and (az.min() < 1 or az.max() > nl or r.min() < 1 or r.max() > width):
        raise ValueError("pixel coordinates outside the raster")
    flat = (az.ravel() - 1) * width + (r.ravel() - 1)

    if rf.is_complex and rf.storage.kind != "f":
        raw = np.memmap(os.fspath(path), dtype=rf.storage, mode="r", shape=(nl, 2 * width)).reshape(-1)
        out = np.empty(flat.size, dtype=rf.dtype)
        for i in range(0, flat.size, chunk):
            idx = flat[i : i + chunk]
            out[i : i + chunk] = raw[2 * idx] + 1j * raw[2 * idx + 1]
        return out.reshape(az.shape)

    raw = np.memmap(os.fspath(path), dtype=rf.dtype, mode="r", shape=(nl, width)).reshape(-1)
    out = np.empty(flat.size, dtype=rf.dtype)
    for i in range(0, flat.size, chunk):
        out[i : i + chunk] = raw[flat[i : i + chunk]]
    return out.reshape(az.shape)
