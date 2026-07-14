from __future__ import annotations

from pathlib import Path
import json
import sys

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_PATH = REPOSITORY_ROOT / "benchmarks" / "reference" / "quantum_v0_1.json"
sys.path.insert(0, str(REPOSITORY_ROOT))

from benchmarks.reference.generate_quantum_v0_1 import collect_reference_metrics


def test_quantum_v0_1_reference_metrics() -> None:
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))

    assert reference["baseline"] == "quantum-v0.1"
    assert reference["backend"] == "numpy"
    assert reference["schema_version"] == 1

    actual_cases = collect_reference_metrics(REPOSITORY_ROOT)
    expected_cases = reference["cases"]
    assert actual_cases.keys() == expected_cases.keys()

    tolerance = reference["tolerance"]
    for case_name, expected_case in expected_cases.items():
        actual_case = actual_cases[case_name]
        assert actual_case["config"] == expected_case["config"]
        assert actual_case["variant"] == expected_case["variant"]
        assert actual_case["metrics"].keys() == expected_case["metrics"].keys()
        for metric_name, expected_value in expected_case["metrics"].items():
            assert actual_case["metrics"][metric_name] == pytest.approx(
                expected_value,
                abs=tolerance["absolute"],
                rel=tolerance["relative"],
            ), f"{case_name}.{metric_name} changed from the quantum v0.1 baseline"
