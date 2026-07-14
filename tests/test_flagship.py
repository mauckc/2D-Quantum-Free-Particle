from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest


jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)

REPOSITORY_ROOT = Path(__file__).parents[1]
REFERENCE_PATH = REPOSITORY_ROOT / "benchmarks" / "reference" / "flagship_m4.json"
sys.path.insert(0, str(REPOSITORY_ROOT))

from benchmarks.reference.generate_flagship_m4 import generate_metrics
from quantum_dynamics_lab.flagship import BASELINE_NAMES, run_flagship
from quantum_dynamics_lab.flagship_config import load_flagship_config


CONFIG_PATH = REPOSITORY_ROOT / "examples" / "robust_flagship.toml"


def test_flagship_scenario_sets_and_required_baselines() -> None:
    config = load_flagship_config(CONFIG_PATH)
    optimization_signatures = {
        scenario.signature for scenario in config.optimization_scenarios
    }
    held_signatures = {scenario.signature for scenario in config.held_out_scenarios}
    required_kinds = {"wavelength", "alignment", "phase_depth", "plane_spacing"}

    assert optimization_signatures.isdisjoint(held_signatures)
    assert required_kinds <= {
        scenario.kind for scenario in config.optimization_scenarios
    }
    assert required_kinds <= {scenario.kind for scenario in config.held_out_scenarios}
    assert BASELINE_NAMES == (
        "no_masks",
        "one_optimized_mask",
        "three_nominal_masks",
        "three_robust_masks",
    )


def test_flagship_reference_is_reproducible_and_reports_all_comparisons() -> None:
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    actual = generate_metrics()

    assert actual["required_baselines"] == list(BASELINE_NAMES)
    assert actual["scenario_sets_disjoint"]
    assert set(actual["baselines"]) == set(BASELINE_NAMES)
    for name in BASELINE_NAMES:
        actual_baseline = actual["baselines"][name]
        expected_baseline = reference["baselines"][name]
        assert len(actual_baseline["held_out"]) == len(actual["held_out_scenarios"])
        for section in ("nominal", "held_out_mean", "held_out_worst"):
            for metric in (
                "mode_overlap",
                "intensity_similarity",
                "mode_power_efficiency",
            ):
                assert actual_baseline[section][metric] == pytest.approx(
                    expected_baseline[section][metric],
                    abs=2e-10,
                    rel=2e-10,
                )
    assert actual["robust_worst_case_overlap_advantage"] == pytest.approx(
        reference["robust_worst_case_overlap_advantage"],
        abs=2e-10,
        rel=2e-10,
    )
    assert actual["robust_worst_case_overlap_advantage"] > reference[
        "declared_limits"
    ]["minimum_robust_worst_case_advantage"]


def test_flagship_writes_comparison_artifacts(tmp_path: Path) -> None:
    result = run_flagship(load_flagship_config(CONFIG_PATH), tmp_path)

    assert result.metrics["robustness_claim_passed"]
    for name in (
        "designs.npz",
        "environment.json",
        "histories.json",
        "metrics.json",
        "report.md",
        "resolved_config.json",
        "scenarios.json",
    ):
        assert (tmp_path / name).exists()
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    for baseline in BASELINE_NAMES:
        assert baseline in report
