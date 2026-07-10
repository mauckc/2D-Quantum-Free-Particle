from __future__ import annotations

from pathlib import Path
import argparse

from .artifacts import save_run, write_metrics
from .config import apply_sweep_variant, load_experiment_config, load_sweep_config
from .experiments import run_experiment
from .render import render_comparison, render_research_narrative, render_run
from .validation import load_validation_config, run_validation


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="quantum-lab")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run one experiment config")
    run_parser.add_argument("config", type=Path)
    run_parser.add_argument("--out", type=Path, required=True)
    run_parser.add_argument("--skip-render", action="store_true")

    render_parser = subparsers.add_parser("render", help="Render an existing run.npz")
    render_parser.add_argument("run", type=Path)
    render_parser.add_argument("--out", type=Path, required=True)

    compare_parser = subparsers.add_parser("compare", help="Run and render a sweep")
    compare_parser.add_argument("config", type=Path)
    compare_parser.add_argument("--out", type=Path, required=True)
    compare_parser.add_argument("--skip-render", action="store_true")

    narrative_parser = subparsers.add_parser(
        "narrative", help="Build a research narrative from generated artifacts"
    )
    narrative_parser.add_argument("--zeno", type=Path, required=True)
    narrative_parser.add_argument("--double-slit", type=Path, required=True)
    narrative_parser.add_argument("--out", type=Path, required=True)

    validate_parser = subparsers.add_parser(
        "validate", help="Run physics validation benchmarks"
    )
    validate_parser.add_argument("config", type=Path)
    validate_parser.add_argument("--out", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "run":
        config = load_experiment_config(args.config)
        result = run_experiment(config)
        run_path = save_run(result, args.out)
        if not args.skip_render:
            render_run(run_path, args.out)
        print(run_path)
        return 0
    if args.command == "render":
        render_run(args.run, args.out)
        print(args.out)
        return 0
    if args.command == "compare":
        sweep = load_sweep_config(args.config)
        base = load_experiment_config(sweep.base_config)
        run_dirs: list[Path] = []
        summary = {"name": sweep.name, "variants": []}
        for variant in sweep.variants:
            config = apply_sweep_variant(base, variant)
            result = run_experiment(config)
            run_dir = args.out / "runs" / variant.name
            run_path = save_run(result, run_dir)
            if not args.skip_render:
                render_run(run_path, run_dir)
            run_dirs.append(run_dir)
            summary["variants"].append(result.metrics)
        if not args.skip_render:
            render_comparison(run_dirs, args.out)
        write_metrics(args.out / "metrics.json", summary)
        print(args.out)
        return 0
    if args.command == "narrative":
        render_research_narrative(args.zeno, args.double_slit, args.out)
        print(args.out)
        return 0
    if args.command == "validate":
        config = load_validation_config(args.config)
        metrics = run_validation(config, args.out)
        print(args.out)
        return 0 if metrics["passed"] else 1
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
