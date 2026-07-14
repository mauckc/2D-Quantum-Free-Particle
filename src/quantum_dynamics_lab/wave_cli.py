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

    optimize_parser = subparsers.add_parser(
        "optimize",
        help="Run constrained phase-mask inverse design",
    )
    optimize_parser.add_argument("config", type=Path)
    optimize_parser.add_argument("--out", type=Path, required=True)
    optimize_parser.add_argument("--resume", type=Path)

    flagship_parser = subparsers.add_parser(
        "flagship",
        help="Run the robust Gaussian-to-HG10 baseline comparison",
    )
    flagship_parser.add_argument("config", type=Path)
    flagship_parser.add_argument("--out", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "propagate":
        config = load_wave_config(args.config)
        run_path = save_wave_run(run_wave_propagation(config), args.out)
        print(run_path)
        return 0
    if args.command == "optimize":
        from .design_config import load_inverse_design_config
        from .optimization import run_inverse_design

        config = load_inverse_design_config(args.config)
        result = run_inverse_design(
            config,
            args.out,
            resume=args.resume,
        )
        print(result.out_dir)
        return 0 if result.status in {"complete", "early_stopped"} else 1
    if args.command == "flagship":
        from .flagship import run_flagship
        from .flagship_config import load_flagship_config

        run_flagship(load_flagship_config(args.config), args.out)
        print(args.out)
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
