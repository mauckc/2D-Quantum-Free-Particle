from __future__ import annotations

import json
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).parents[1]
REFERENCE_PATH = REPOSITORY_ROOT / "benchmarks" / "reference" / "optics_m1.json"
sys.path.insert(0, str(REPOSITORY_ROOT))

from benchmarks.reference.generate_optics_m1 import generate_metrics


def test_optics_m1_reference_metrics_meet_declared_limits() -> None:
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    regenerated = generate_metrics()

    assert reference["schema_version"] == 1
    assert regenerated["configuration"] == reference["configuration"]
    metrics = regenerated["metrics"]
    limits = reference["declared_limits"]
    assert metrics["gaussian_relative_l2_fine"] < limits[
        "gaussian_relative_l2_fine_max"
    ]
    assert metrics["gaussian_convergence_ratio"] < limits[
        "gaussian_convergence_ratio_max"
    ]
    assert metrics["focal_waist_relative_error"] < limits[
        "focal_waist_relative_error_max"
    ]
    assert metrics["fresnel_power_relative_drift"] < limits[
        "power_relative_drift_max"
    ]
    assert metrics["angular_power_relative_drift"] < limits[
        "power_relative_drift_max"
    ]
    assert metrics["fresnel_reverse_max_error"] < limits["reverse_max_error_max"]
    assert metrics["angular_reverse_max_error"] < limits["reverse_max_error_max"]
    assert metrics["cross_model_complex_relative_l2"] < limits[
        "cross_model_relative_l2_max"
    ]
    assert metrics["cross_model_intensity_relative_l2"] < limits[
        "cross_model_relative_l2_max"
    ]
    assert metrics["angular_retained_spectral_power"] > limits[
        "minimum_retained_spectral_power"
    ]
