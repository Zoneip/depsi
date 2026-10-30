"""Command line entry point: ``depsipy run <param_file.txt>``.

The CLI is deliberately thin: it parses the requested processing groups and
delegates to :func:`depsipy.pipeline.run_pipeline`, which writes every DePSI
product into the work directory.
"""

from __future__ import annotations

import argparse
import json
import sys

from . import DEPSI_COMPAT, __version__
from .pipeline import parse_groups, run_pipeline

__all__ = ["main", "build_parser"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="depsipy",
        description=(
            "DePSI-compatible Persistent Scatterer InSAR time-series processor "
            f"(models DePSI {DEPSI_COMPAT})."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the processing chain for a DePSI parameter file")
    run.add_argument("param_file", help="path to param_file.txt")
    run.add_argument(
        "--workdir",
        default=None,
        help="directory for the products (default: the parameter file's directory)",
    )
    run.add_argument(
        "--groups",
        default=None,
        help="processing groups to run, e.g. '3,4' or '[1:8]' (default: 1..8)",
    )
    run.add_argument(
        "--json",
        action="store_true",
        help="print the run summary as JSON instead of key/value lines",
    )
    run.add_argument("--quiet", action="store_true", help="suppress the progress lines")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    def progress(stage: str, elapsed: float) -> None:
        if not args.quiet:
            print(f"[{elapsed:8.2f}s] {stage}", file=sys.stderr, flush=True)

    groups = parse_groups(args.groups)
    summary = run_pipeline(
        args.param_file,
        workdir=args.workdir,
        groups=groups,
        progress=progress,
    )
    if args.json:
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        for key, value in summary.items():
            print(f"{key}: {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
