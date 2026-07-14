from __future__ import annotations

import json
from pathlib import Path

import pytest


jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)

from quantum_dynamics_lab.flagship import BASELINE_NAMES, run_flagship
from quantum_dynamics_lab.flagship_config import load_flagship_config
from quantum_dynamics_lab.flagship_validation import validate_flagship


REPOSITORY_ROOT = Path(__file__).parents[1]
CONFIG_PATH = REPOSITORY_ROOT / "examples" / "robust_flagship.toml"
REFERENCE_PATH = (
    REPOSITORY_ROOT
    / "benchmarks"
    / "reference"
    / "flagship_validation_m4.json"
)


def test_flagship_model_transfer_and_convergence_reference(tmp_path: Path) -> None:
    config = load_flagship_config(CONFIG_PATH)
    flagship = run_flagship(config, tmp_path)
    actual = validate_flagship(config, flagship, tmp_path).metrics
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    tolerance = actual["declared_limits"]

    assert actual["trusted"]
    assert all(actual["checks"].values())
    assert actual["checks"] == reference["checks"]
    for baseline in BASELINE_NAMES:
        comparisons = actual["model_transfer"][baseline]
        assert len(comparisons) == 1 + len(config.held_out_scenarios)
        assert all("fresnel" in item and "angular_spectrum" in item for item in comparisons)
    for key, expected in reference["model_transfer_summary"].items():
        assert actual["model_transfer_summary"][key] == pytest.approx(
            expected,
            abs=tolerance["reproduction_absolute_tolerance"],
            rel=tolerance["reproduction_relative_tolerance"],
        )
    convergence = actual["convergence"]
    for key, expected in reference["convergence_summary"].items():
        assert convergence[key] == pytest.approx(
            expected,
            abs=tolerance["reproduction_absolute_tolerance"],
            rel=tolerance["reproduction_relative_tolerance"],
        )
    for name, key in (
        ("grid_overlap", "grid"),
        ("padding_overlap", "padding"),
        ("aperture_overlap", "aperture"),
        ("step_overlap", "step"),
    ):
        assert [item["mode_overlap"] for item in convergence[key]] == pytest.approx(
            reference["series"][name],
            abs=tolerance["reproduction_absolute_tolerance"],
            rel=tolerance["reproduction_relative_tolerance"],
        )
    assert not actual["maxwell_fdtd"]["performed"]
    assert not actual["maxwell_fdtd"]["required"]
    assert (tmp_path / "validation_metrics.json").exists()
    assert (tmp_path / "model_transfer_report.md").exists()
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "Independent validation" in report
    assert "Trusted: **True**" in report
