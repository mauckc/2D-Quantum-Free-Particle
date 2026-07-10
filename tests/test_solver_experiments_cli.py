from __future__ import annotations

import numpy as np

from quantum_dynamics_lab.cli import main
from quantum_dynamics_lab.config import (
    BoundaryConfig,
    ExperimentConfig,
    GridConfig,
    MeasurementConfig,
    PotentialConfig,
    SolverConfig,
    WavePacketConfig,
    load_sweep_config,
)
from quantum_dynamics_lab.experiments import run_experiment
from quantum_dynamics_lab.grid import make_grid, norm, normalize
from quantum_dynamics_lab.render import render_run
from quantum_dynamics_lab.solver import SplitStepSolver


def test_free_packet_norm_drift_stays_small() -> None:
    result = run_experiment(
        ExperimentConfig(
            name="free_norm",
            experiment="free_packet",
            grid=GridConfig(n=32, length=8.0),
            wave_packet=WavePacketConfig(
                center=(-1.5, 0.0), sigma=0.5, momentum=(2.0, 0.5)
            ),
            potential=PotentialConfig(kind="free"),
            solver=SolverConfig(dt=0.002, tf=0.08, frame_interval=0.02),
        )
    )

    drift = np.max(np.abs(result.arrays["norm"] - 1.0))
    assert drift < 1e-8


def test_plane_wave_accumulates_expected_free_phase() -> None:
    grid = make_grid(GridConfig(n=32, length=8.0))
    mode = 2
    kx = 2.0 * np.pi * mode / grid.length
    psi0 = normalize(np.exp(1j * kx * grid.X).astype(np.complex128), grid)
    solver_config = SolverConfig(dt=0.01, tf=0.01, frame_interval=0.01)
    solver = SplitStepSolver(grid, np.zeros((grid.n, grid.n)), solver_config)

    evolved = solver.step(psi0)
    expected = psi0 * np.exp(-1j * kx * kx * solver_config.dt / 2.0)

    assert np.max(np.abs(evolved - expected)) < 1e-12
    assert np.isclose(norm(evolved, grid), 1.0)


def test_barrier_zeno_fast_measurement_reduces_detection() -> None:
    common = dict(
        experiment="barrier_zeno",
        grid=GridConfig(n=40, length=8.0),
        wave_packet=WavePacketConfig(center=(-2.4, 0.0), sigma=0.45, momentum=(4.5, 0.0)),
        potential=PotentialConfig(
            kind="barrier", height=35.0, barrier_x=0.0, barrier_width=0.25
        ),
        solver=SolverConfig(dt=0.002, tf=0.9, frame_interval=0.03),
    )
    slow = run_experiment(
        ExperimentConfig(
            name="slow",
            measurement=MeasurementConfig(kind="zeno", interval=0.24),
            **common,
        )
    )
    fast = run_experiment(
        ExperimentConfig(
            name="fast",
            measurement=MeasurementConfig(kind="zeno", interval=0.06),
            **common,
        )
    )

    assert (
        fast.metrics["final_detected_probability"]
        < slow.metrics["final_detected_probability"]
    )
    assert "final_unconditional_transmission_probability" in fast.metrics


def test_absorbing_boundary_tracks_probability_loss() -> None:
    result = run_experiment(
        ExperimentConfig(
            name="absorber",
            experiment="free_packet",
            grid=GridConfig(n=32, length=8.0),
            wave_packet=WavePacketConfig(center=(2.6, 0.0), sigma=0.45, momentum=(5.0, 0.0)),
            potential=PotentialConfig(kind="free"),
            solver=SolverConfig(dt=0.002, tf=0.25, frame_interval=0.05),
            boundary=BoundaryConfig(kind="absorbing", width=1.2, strength=5.0),
        )
    )

    assert result.metrics["final_absorbed_probability"] > 0.0
    assert result.metrics["final_survival_weight"] < 1.0


def test_double_slit_which_path_reduces_interference_contrast() -> None:
    result = run_experiment(
        ExperimentConfig(
            name="double_slit_smoke",
            experiment="double_slit",
            grid=GridConfig(n=40, length=8.0),
            wave_packet=WavePacketConfig(center=(-2.4, 0.0), sigma=0.45, momentum=(5.0, 0.0)),
            potential=PotentialConfig(
                kind="double_slit",
                height=100.0,
                barrier_x=0.0,
                barrier_width=0.25,
                slit_width=0.45,
                slit_separation=1.1,
            ),
            solver=SolverConfig(dt=0.002, tf=0.9, frame_interval=0.03),
            measurement=MeasurementConfig(kind="which_path"),
        )
    )

    assert (
        result.metrics["which_path_interference_contrast"]
        < result.metrics["coherent_interference_contrast"]
    )


def test_cli_run_and_render_smoke(tmp_path) -> None:
    out = tmp_path / "free_packet"
    assert (
        main(
            [
                "run",
                "examples/free_packet.toml",
                "--out",
                str(out),
                "--skip-render",
            ]
        )
        == 0
    )
    assert (out / "run.npz").exists()
    assert main(["render", str(out / "run.npz"), "--out", str(out / "report")]) == 0
    assert (out / "report" / "summary.png").exists()
    assert (out / "report" / "potential.png").exists()
    assert (out / "report" / "report.md").exists()


def test_research_sweep_scan_expands_variants() -> None:
    sweep = load_sweep_config("examples/zeno_research_sweep.toml")

    assert len(sweep.variants) == 12
    assert sweep.variants[0].potential_height == 25.0
    assert sweep.variants[0].measurement_interval is None


def test_compare_generates_polished_report(tmp_path) -> None:
    base = tmp_path / "base.toml"
    base.write_text(
        """
name = "small_zeno"
experiment = "barrier_zeno"

[grid]
n = 24
length = 8.0

[wave_packet]
center = [-2.0, 0.0]
sigma = 0.45
momentum = [4.0, 0.0]

[potential]
kind = "barrier"
height = 30.0
barrier_x = 0.0
barrier_width = 0.25

[solver]
dt = 0.004
tf = 0.12
frame_interval = 0.04

[boundary]
kind = "absorbing"
width = 1.0
strength = 4.0

[measurement]
kind = "zeno"
interval = 0.08
detector_buffer = 0.25

[output]
render_video = false
""",
        encoding="utf-8",
    )
    sweep = tmp_path / "sweep.toml"
    sweep.write_text(
        """
name = "small_sweep"
base_config = "base.toml"

[scan]
potential_heights = [30.0]
measurement_intervals = [0.0, 0.08]
""",
        encoding="utf-8",
    )
    out = tmp_path / "comparison"

    assert main(["compare", str(sweep), "--out", str(out)]) == 0
    assert (out / "comparison.csv").exists()
    assert (out / "comparison.png").exists()
    assert (out / "zeno_transmission_heatmap.png").exists()
    assert (out / "report.md").exists()
