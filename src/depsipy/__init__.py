"""depsipy -- a fast, feature-complete reimplementation of the DePSI
Persistent Scatterer InSAR time-series toolbox (TU Delft, v1.7.8).

The package mirrors the module structure and file formats of the MATLAB
toolbox so that existing projects can be processed unchanged, while the
numerical core is written with numpy/scipy (+ optional numba) kernels.

Public entry points
-------------------
``depsipy.io``      raw binary / parameter-file I/O (drop-in for freadbk/readinput)
``depsipy.core``    functional models, LAMBDA integer least squares, statistics
``depsipy.steps``   the processing chain, one module per DePSI processing group
``depsipy.cli``     ``depsipy run param_file.txt``
"""

from __future__ import annotations

__version__ = "0.1.0"

# DePSI tool version this implementation is modelled on.
DEPSI_COMPAT = "1.7.8.8"

__all__ = ["__version__", "DEPSI_COMPAT"]
