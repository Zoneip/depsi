"""I/O subpackage: DePSI raw rasters, parameter files and project state."""

from __future__ import annotations

from .doris import (
    DorisMetadata,
    DorisRaster,
    MFFHeader,
    doris_format,
    read_doris_metadata,
    read_mff_header,
    read_res,
)
from .raw import (
    RawFile,
    RawFormat,
    file_size,
    infer_shape,
    iter_blocks,
    parse_format,
    read_pixels,
    read_raw,
    write_raw,
)

__all__ = [
    "DorisMetadata",
    "DorisRaster",
    "MFFHeader",
    "doris_format",
    "read_doris_metadata",
    "read_mff_header",
    "read_res",
    "RawFile",
    "RawFormat",
    "file_size",
    "infer_shape",
    "iter_blocks",
    "parse_format",
    "read_pixels",
    "read_raw",
    "write_raw",
]
