from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest


jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)

REPOSITORY_ROOT = Path(__file__).parents[1]
REFERENCE_PATH = REPOSITORY_ROOT / "benchmarks" / "reference" / "jax_m2.json"
sys.path.insert(0, str(REPOSITORY_ROOT))

from benchmarks.reference.generate_jax_m2 import generate_metrics


def test_jax_objective_parity_gradients_and_compiled_update() -> None:
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    evidence = generate_metrics()
    limits = reference["declared_limits"]
    metrics = evidence["metrics"]

    assert evidence["seed"] == reference["seed"]
    assert evidence["configuration"] == reference["configuration"]
    assert metrics["objective_value_max_abs_error"] < limits[
        "objective_value_max_abs_error"
    ]
    for name in ("mode_gradient", "intensity_gradient"):
        assert abs(metrics[name]["autodiff_directional_derivative"]) > 1e-12
        assert metrics[name]["absolute_error"] < limits["gradient_absolute_error"]
        assert metrics[name]["relative_error"] < limits["gradient_relative_error"]
    assert metrics["compiled_update_all_finite"]
    assert metrics["compiled_update_gradient_norm"] > 0.0
    assert jax.default_backend() == "cpu"
