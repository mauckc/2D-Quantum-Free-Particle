from __future__ import annotations

import json

from quantum_dynamics_lab.cli import main
from quantum_dynamics_lab.validation import load_validation_config, run_validation


def test_validation_config_loads_example() -> None:
    config = load_validation_config("examples/validation_suite.toml")

    assert config.name == "physics_validation_suite"
    assert len(config.norm_cases) == 2
    assert config.barrier_case is not None
    assert config.zeno_case is not None
    assert config.backend_case is not None


def test_validation_runner_writes_artifacts(tmp_path) -> None:
    config_path = tmp_path / "validation.toml"
    config_path.write_text(
        """
name = "small_validation"

[[norm_cases]]
name = "norm_small"
n = 24
dt = 0.004
tf = 0.04
tolerance = 1e-8

[[dispersion_cases]]
name = "dispersion_small"
n = 32
dt = 0.004
tf = 0.04
sigma = 0.7
tolerance = 0.5

[barrier_case]
name = "barrier_small"
heights = [20.0, 40.0]
tolerance = 1e-6

[zeno_case]
name = "zeno_small"
intervals = [0.0, 0.06, 0.12]
tolerance = 1e-6

[backend_case]
name = "backend_small"
backends = ["numpy", "pyfftw", "cupy"]
tolerance = 1e-10
""",
        encoding="utf-8",
    )
    out = tmp_path / "validation"

    assert main(["validate", str(config_path), "--out", str(out)]) == 0

    metrics = json.loads((out / "validation_metrics.json").read_text(encoding="utf-8"))
    assert metrics["passed"] is True
    assert {result["kind"] for result in metrics["results"]} == {
        "norm_conservation",
        "free_gaussian_dispersion",
        "barrier_transmission_monotonicity",
        "zeno_no_click_transmission",
        "backend_parity",
    }
    assert (out / "validation_report.md").exists()
    assert (out / "validation_report.html").exists()
    assert (out / "validation_summary.png").exists()
    assert (out / "norm_convergence.png").exists()
    assert (out / "dispersion_error.png").exists()
    assert (out / "barrier_transmission.png").exists()
    assert (out / "zeno_transmission.png").exists()


def test_run_validation_returns_failed_status_for_tight_tolerance(tmp_path) -> None:
    config_path = tmp_path / "validation.toml"
    config_path.write_text(
        """
name = "expected_failure"

[[dispersion_cases]]
name = "too_tight"
n = 24
dt = 0.004
tf = 0.08
sigma = 0.7
tolerance = 0.0
""",
        encoding="utf-8",
    )
    config = load_validation_config(config_path)
    metrics = run_validation(config, tmp_path / "out")

    assert metrics["passed"] is False
