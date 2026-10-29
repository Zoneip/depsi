"""DePSI parameter-file parsing (drop-in for ``readinput.m``).

The MATLAB implementation tokenises the whole file with ``textread(f,'%s')``
and then takes the token *two positions after* the parameter name (the ``=``
sits in between), finally interpreting the token with ``str2num`` (which is
really ``eval``).  That makes ``[1:9]``, ``[30,100]``, ``2*pi`` and ``[]`` all
valid values, and it also means a value spread over several whitespace
separated tokens confuses it.

This parser is line/assignment oriented, which is strictly more forgiving:
values may span lines, trailing ``;`` and ``%`` comments are ignored, MATLAB
range literals and simple arithmetic are evaluated with a whitelist
(no ``eval``).
"""

from __future__ import annotations

import ast
import math
import os
import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np

__all__ = [
    "ParamError",
    "parse_param_file",
    "parse_value",
    "Params",
    "load_params",
]

_ASSIGN_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.*)$")
_IDENT_RE = re.compile(r"^[A-Za-z_][A-Za-z_0-9]*$")


class ParamError(ValueError):
    """Raised for malformed or missing parameters."""


# ---------------------------------------------------------------------------
# MATLAB expression evaluation (safe subset)
# ---------------------------------------------------------------------------

_ALLOWED_NODES = (
    ast.Expression,
    ast.BinOp,
    ast.UnaryOp,
    ast.Constant,
    ast.List,
    ast.Tuple,
    ast.Call,
    ast.Name,
    ast.Load,
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.Pow,
    ast.USub,
    ast.UAdd,
)

_MATH_NAMES = {
    "pi": math.pi,
    "e": math.e,
    "inf": math.inf,
    "Inf": math.inf,
    "nan": math.nan,
    "NaN": math.nan,
    "sqrt": math.sqrt,
    "abs": abs,
    "exp": math.exp,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "log": math.log,
    "log10": math.log10,
}


def _strip_comments(text: str) -> str:
    """Drop a trailing ``%`` comment, respecting single quotes."""
    out = []
    in_str = False
    for ch in text:
        if ch == "'":
            in_str = not in_str
        elif ch == "%" and not in_str:
            break
        out.append(ch)
    return "".join(out)


def parse_value(text: str) -> Any:
    """Interpret the right-hand side of a parameter assignment.

    Returns ``float``/``int`` scalars, ``numpy`` arrays for ``[...]``, or a
    ``str`` for quoted values (quotes removed, as ``readinput.m`` does).
    """
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    s = _strip_comments(s).strip()
    if s == "":
        return ""

    # MATLAB char literal: 'text'  (doubled quotes are escapes)
    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
        return s[1:-1].replace("''", "'")

    # MATLAB numeric literal array: [a, b, c]
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        return _parse_array_literal(inner)

    val = _eval_expr(s)
    return _normalise_scalar(val)


def _normalise_scalar(val: Any) -> Any:
    if isinstance(val, (list, tuple, np.ndarray)):
        return np.asarray(val)
    if isinstance(val, complex):
        return val
    if isinstance(val, bool):
        return val
    if isinstance(val, float) and val.is_integer() and abs(val) < 2**53:
        # DePSI never distinguishes 1 from 1.0; ints behave better as indices
        return int(val) if not isinstance(val, np.floating) else val
    return val


def _parse_array_literal(inner: str) -> Any:
    inner = inner.strip()
    if inner == "":
        return np.array([], dtype=float)

    # ranges: a:b or a:step:b  (possibly as one element of a larger list)
    parts = _split_top_level(inner)
    if len(parts) == 1 and parts[0].count(":") in (1, 2):
        return _parse_range(parts[0])

    values: list[Any] = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if part.count(":") in (1, 2) and not re.search(r"[A-Za-z_]", part):
            values.extend(np.atleast_1d(_parse_range(part)).tolist())
        else:
            values.append(_eval_expr(part))
    arr = np.asarray(values)
    if arr.dtype.kind in "fc" and np.all(np.isreal(arr)):
        arr = arr.real
    return arr


def _split_top_level(text: str) -> list[str]:
    """Split on commas / whitespace at bracket depth 0, keeping ranges intact."""
    parts: list[str] = []
    depth = 0
    cur = ""
    for ch in text:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if depth == 0 and (ch == "," or ch in " \t\n"):
            if cur.strip():
                parts.append(cur.strip())
            cur = "" if ch != "," else ""
            continue
        cur += ch
    if cur.strip():
        parts.append(cur.strip())
    return parts


def _parse_range(text: str) -> np.ndarray:
    tokens = [t.strip() for t in text.split(":")]
    if len(tokens) == 2:
        start, stop = (_eval_expr(tokens[0]), _eval_expr(tokens[1]))
        step = 1.0
    elif len(tokens) == 3:
        start, step, stop = (_eval_expr(t) for t in tokens)
    else:
        raise ParamError(f"cannot parse range {text!r}")
    if step == 0:
        raise ParamError(f"zero step in range {text!r}")
    count = int(math.floor((stop - start) / step + 1e-9)) + 1
    count = max(count, 0)
    arr = start + step * np.arange(count)
    if float(arr[-1]) > stop if len(arr) else False:
        arr = arr[:-1]
    if np.all(np.isreal(arr)) and np.all(arr == np.round(arr)):
        return arr.astype(np.int64)
    return arr.astype(float)


def _eval_expr(text: str) -> Any:
    text = text.strip()
    if not text:
        return np.nan
    # fast path: plain number
    try:
        return float(text) if any(c in text for c in ".eE") else int(text)
    except ValueError:
        pass
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as exc:  # pragma: no cover - defensive
        raise ParamError(f"cannot parse value {text!r}: {exc}") from exc
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise ParamError(f"unsupported expression {text!r}")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in _MATH_NAMES:
                raise ParamError(f"unsupported function in {text!r}")
    try:
        return eval(compile(tree, "<param>", "eval"), {"__builtins__": {}}, dict(_MATH_NAMES))
    except Exception as exc:
        raise ParamError(f"cannot evaluate {text!r}: {exc}") from exc


# ---------------------------------------------------------------------------
# file parsing
# ---------------------------------------------------------------------------


def parse_param_file(path: str | os.PathLike) -> dict[str, Any]:
    """Parse a DePSI ``param_file.txt`` into an ordered ``{name: value}`` dict."""
    with open(os.fspath(path), "r", encoding="utf-8", errors="replace") as fh:
        raw_lines = fh.readlines()

    result: dict[str, Any] = {}
    i = 0
    n = len(raw_lines)
    while i < n:
        line = _strip_comments(raw_lines[i])
        m = _ASSIGN_RE.match(line)
        if not m:
            i += 1
            continue
        name, rhs = m.group(1), m.group(2)
        # accumulate continuation lines while brackets are unbalanced
        joined = rhs.rstrip("\n")
        j = i
        while _bracket_balance(joined) > 0 and j + 1 < n:
            j += 1
            joined += " " + _strip_comments(raw_lines[j]).strip()
        if name in result:
            raise ParamError(f"parameter {name!r} appears more than once in {path}")
        try:
            result[name] = parse_value(joined)
        except ParamError as exc:
            raise ParamError(f"{path}: parameter {name!r}: {exc}") from exc
        i = j + 1
    return result


def _bracket_balance(text: str) -> int:
    depth = 0
    in_str = False
    for ch in text:
        if ch == "'":
            in_str = not in_str
        elif not in_str:
            if ch in "[(":
                depth += 1
            elif ch in "])":
                depth -= 1
    return depth


# ---------------------------------------------------------------------------
# typed accessor
# ---------------------------------------------------------------------------


@dataclass
class Params:
    """Typed, validated view over a DePSI parameter dictionary.

    Defaults follow ``ps_readinput_parameters.m``; a ``ParamError`` is raised
    for unknown defaults only when the parameter is actually used.
    """

    values: dict[str, Any] = field(default_factory=dict)
    path: str | None = None

    # -- access ------------------------------------------------------------
    def get(self, name: str, default: Any = ...):
        if name in self.values:
            return self.values[name]
        if default is ...:
            raise ParamError(f"required parameter {name!r} not found in parameter file")
        return default

    def has(self, name: str) -> bool:
        return name in self.values

    def __contains__(self, name: str) -> bool:
        return name in self.values

    def __getitem__(self, name: str) -> Any:
        return self.get(name)

    # -- typed helpers -----------------------------------------------------
    def int(self, name: str, default=...) -> int:
        v = self.get(name, default)
        return int(np.asarray(v).ravel()[0])

    def float(self, name: str, default=...) -> float:
        v = self.get(name, default)
        return float(np.asarray(v).ravel()[0])

    def str(self, name: str, default=...) -> str:
        v = self.get(name, default)
        if isinstance(v, str):
            return v
        raise ParamError(f"parameter {name!r} should be a string, got {v!r}")

    def ints(self, name: str, default=...) -> np.ndarray:
        v = self.get(name, default)
        return np.atleast_1d(np.asarray(v)).astype(np.int64)

    def floats(self, name: str, default=...) -> np.ndarray:
        v = self.get(name, default)
        return np.atleast_1d(np.asarray(v)).astype(float)

    def is_empty(self, name: str) -> bool:
        v = self.values.get(name, None)
        if v is None or (isinstance(v, str) and v == ""):
            return True
        arr = np.asarray(v)
        return arr.size == 0

    def as_dict(self) -> dict[str, Any]:
        return dict(self.values)


def load_params(path: str | os.PathLike) -> Params:
    return Params(values=parse_param_file(path), path=os.fspath(path))
