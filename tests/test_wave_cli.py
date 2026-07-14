from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from quantum_dynamics_lab.cli import main as quantum_main
from quantum_dynamics_lab.wave_cli import main as wave_main
from quantum_dynamics_lab.wave_config import load_wave_config
from quantum_dynamics_lab.wave_experiments import run_wave_propagation


REPOSITORY_ROOT = Path(__file__).parents[1]


def test_committed_gaussian_lens_config_reaches_predicted_focus(tmp_path: Path) -> None:
    config_path = REPOSITORY_ROOT / "examples" / "gaussian_lens.toml"
    config = load_wave_config(config_path)
    run = run_wave_propagation(config)
    expected_waist = (
        config.grid.wavelength_m
        * config.steps[0].focal_length_m
        / (math.pi * config.source.waist_m)
    )

    assert run.metrics["output_second_moment_waist_m"][1] == pytest.approx(
        expected_waist,
        rel=0.01,
    )
    assert wave_main(["propagate", str(config_path), "--out", str(tmp_path)]) == 0
    assert (tmp_path / "wave_run.npz").exists()
    metrics = json.loads((tmp_path / "metrics.json").read_text(encoding="utf-8"))
    resolved = json.loads(
        (tmp_path / "resolved_config.json").read_text(encoding="utf-8")
    )
    assert metrics["name"] == "gaussian_lens_focus"
    assert metrics["input_power"] == pytest.approx(1.0, abs=1e-13)
    assert metrics["final_power"] == pytest.approx(1.0, abs=1e-13)
    assert [stage["kind"] for stage in metrics["stages"]] == ["lens", "propagate"]
    assert resolved["grid"]["wavelength_m"] == config.grid.wavelength_m
    with np.load(tmp_path / "wave_run.npz") as archive:
        assert archive["input_field"].shape == config.grid.shape
        assert archive["final_field"].shape == config.grid.shape
        metadata = json.loads(str(archive["metadata"]))
    assert metadata["git_commit"]
    assert "git_dirty" in metadata


def test_phase_mask_blas_and_arbitrary_array_source(tmp_path: Path) -> None:
    source = np.ones((32, 48), dtype=np.complex128)
    np.save(tmp_path / "source.npy", source)
    config_path = tmp_path / "wave.toml"
    config_path.write_text(
        """
name = "array_blas"
[grid]
shape = [32, 48]
extent_m = [0.003, 0.004]
wavelength_m = 6.33e-7
[source]
kind = "array"
path = "source.npy"
[[steps]]
kind = "phase_mask"
phase_profile = "sinusoidal"
phase_amplitude_rad = 0.5
phase_period_m = 0.0008
axis = "x"
[[steps]]
kind = "circular_aperture"
radius_m = 0.001
[[steps]]
kind = "propagate"
model = "blas"
distance_m = 0.02
""".strip(),
        encoding="utf-8",
    )

    run = run_wave_propagation(load_wave_config(config_path))

    assert run.input_field.shape == (32, 48)
    assert run.metrics["stages"][2]["retained_spectral_power"] <= 1.0
    assert run.metrics["final_power"] <= run.metrics["input_power"]


def test_invalid_step_and_both_cli_help_surfaces(tmp_path: Path) -> None:
    config_path = tmp_path / "invalid.toml"
    config_path.write_text(
        "name='invalid'\n[[steps]]\nkind='propagate'\nmodel='unknown'\ndistance_m=1.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unsupported propagation model"):
        run_wave_propagation(load_wave_config(config_path))
    with pytest.raises(SystemExit) as wave_help:
        wave_main(["--help"])
    with pytest.raises(SystemExit) as quantum_help:
        quantum_main(["--help"])
    assert wave_help.value.code == 0
    assert quantum_help.value.code == 0
