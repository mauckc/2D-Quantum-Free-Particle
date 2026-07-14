"""Command-line interface for scalar wave propagation and inverse design."""

from __future__ import annotations

import argparse
from pathlib import Path

from .wave_artifacts import save_wave_run
from .wave_config import load_wave_config
from .wave_experiments import run_wave_propagation


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="wave-lab",
        description="Differentiable 2D wave propagation and inverse-design lab",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    propagate_parser = subparsers.add_parser(
        "propagate",
        help="Run a scalar-optics forward config",
    )
    propagate_parser.add_argument("config", type=Path)
    propagate_parser.add_argument("--out", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "propagate":
        config = load_wave_config(args.config)
        run_path = save_wave_run(run_wave_propagation(config), args.out)
        print(run_path)
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
