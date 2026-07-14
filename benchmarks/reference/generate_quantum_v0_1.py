from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any
import argparse
import json
import subprocess

import numpy as np

from quantum_dynamics_lab.config import load_experiment_config
from quantum_dynamics_lab.experiments import ExperimentResult, run_experiment
from quantum_dynamics_lab.grid import integrate_probability


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_PATH = Path(__file__).with_name("quantum_v0_1.json")

CASE_SPECS = (
    ("free_packet", "examples/free_packet.toml", None),
    ("barrier_unmeasured", "examples/barrier_zeno.toml", "unmeasured"),
    ("barrier_zeno", "examples/barrier_zeno.toml", None),
    ("double_slit", "examples/double_slit.toml", None),
)


def collect_reference_metrics(repository_root: Path = REPOSITORY_ROOT) -> dict[str, Any]:
    cases: dict[str, Any] = {}
    for name, config_path, variant in CASE_SPECS:
        config = load_experiment_config(repository_root / config_path)
        config = replace(config, solver=replace(config.solver, backend="numpy"))
        if variant == "unmeasured":
            config = replace(
                config,
                name=name,
                measurement=replace(
                    config.measurement,
                    kind="none",
                    interval=None,
                    time=None,
                ),
            )
        result = run_experiment(config)
        cases[name] = {
            "config": config_path,
            "variant": variant or "committed",
            "metrics": _case_metrics(name, result),
        }
    return cases


def build_reference_document(repository_root: Path = REPOSITORY_ROOT) -> dict[str, Any]:
    return {
        "backend": "numpy",
        "baseline": "quantum-v0.1",
        "cases": collect_reference_metrics(repository_root),
        "milestone": "M0 - Pivot foundation",
        "package_version": "0.1.0",
        "schema_version": 1,
        "scientific_interpretation": (
            "These values freeze numerical behavior. They do not independently "
            "validate the no-click Zeno or global min/max double-slit diagnostics."
        ),
        "source_commit": _source_commit(repository_root),
        "tolerance": {
            "absolute": 1e-10,
            "relative": 1e-9,
            "justification": (
                "The cases use deterministic NumPy complex128 FFT propagation. "
                "The tolerance permits platform-level floating-point roundoff "
                "while remaining sensitive to solver and experiment changes."
            ),
        },
    }


def _case_metrics(name: str, result: ExperimentResult) -> dict[str, float]:
    metrics = {
        "final_norm": float(result.metrics["final_norm"]),
        "max_norm_drift": float(np.max(np.abs(result.arrays["norm"] - 1.0))),
    }
    if name == "free_packet":
        metrics.update(_free_packet_moments(result))
    elif name.startswith("barrier_"):
        for key in (
            "final_absorbed_probability",
            "final_detected_probability",
            "final_survival_weight",
            "final_transmission_probability",
            "final_unconditional_transmission_probability",
        ):
            metrics[key] = float(result.metrics[key])
        metrics["probability_balance"] = (
            metrics["final_absorbed_probability"]
            + metrics["final_detected_probability"]
            + metrics["final_survival_weight"]
        )
    elif name == "double_slit":
        for key in (
            "coherent_interference_contrast",
            "final_which_path_norm",
            "which_path_interference_contrast",
        ):
            metrics[key] = float(result.metrics[key])
        coherent_profile = result.arrays["coherent_screen_profile"]
        which_path_profile = result.arrays["which_path_screen_profile"]
        metrics["coherent_screen_profile_integral"] = float(
            np.sum(coherent_profile) * result.grid.dy
        )
        metrics["which_path_screen_profile_integral"] = float(
            np.sum(which_path_profile) * result.grid.dy
        )
    return metrics


def _free_packet_moments(result: ExperimentResult) -> dict[str, float]:
    probability = np.abs(result.arrays["psi_frames"][-1]) ** 2
    mean_x = integrate_probability(probability * result.grid.X, result.grid)
    mean_y = integrate_probability(probability * result.grid.Y, result.grid)
    return {
        "final_x_mean": mean_x,
        "final_x_variance": integrate_probability(
            probability * (result.grid.X - mean_x) ** 2, result.grid
        ),
        "final_y_mean": mean_y,
        "final_y_variance": integrate_probability(
            probability * (result.grid.Y - mean_y) ** 2, result.grid
        ),
    }


def _source_commit(repository_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the quantum v0.1 metrics")
    parser.add_argument("--output", type=Path, default=REFERENCE_PATH)
    args = parser.parse_args()
    document = build_reference_document()
    args.output.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
