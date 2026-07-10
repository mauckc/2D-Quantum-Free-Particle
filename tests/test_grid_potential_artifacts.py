from __future__ import annotations

import json

import numpy as np

from quantum_dynamics_lab.artifacts import load_run, save_run
from quantum_dynamics_lab.config import (
    ExperimentConfig,
    GridConfig,
    PotentialConfig,
    SolverConfig,
    WavePacketConfig,
)
from quantum_dynamics_lab.experiments import run_experiment
from quantum_dynamics_lab.grid import make_grid, norm
from quantum_dynamics_lab.potentials import (
    build_potential,
    lower_slit_mask,
    transmission_mask,
    upper_slit_mask,
)
from quantum_dynamics_lab.wavepackets import gaussian_packet


def test_grid_spacing_and_fft_wavenumbers() -> None:
    grid = make_grid(GridConfig(n=16, length=8.0))

    assert grid.dx == 0.5
    assert grid.x[0] == -4.0
    assert np.isclose(grid.kx[0, 1], 2.0 * np.pi / 8.0)


def test_wave_packet_is_normalized() -> None:
    grid = make_grid(GridConfig(n=32, length=8.0))
    psi = gaussian_packet(
        WavePacketConfig(center=(-1.0, 0.25), sigma=0.45, momentum=(2.0, -1.0)),
        grid,
    )

    assert np.isclose(norm(psi, grid), 1.0)


def test_potential_masks_describe_barrier_and_slits() -> None:
    grid = make_grid(GridConfig(n=40, length=10.0))
    config = PotentialConfig(
        kind="double_slit",
        height=100.0,
        barrier_x=0.0,
        barrier_width=0.5,
        slit_width=0.5,
        slit_separation=1.5,
    )

    potential = build_potential(config, grid)
    upper = upper_slit_mask(config, grid)
    lower = lower_slit_mask(config, grid)
    transmitted = transmission_mask(config, grid)

    assert potential.shape == (40, 40)
    assert np.any(potential == 100.0)
    assert np.any(upper)
    assert np.any(lower)
    assert not np.any(upper & lower)
    assert np.any(transmitted)


def test_npz_round_trip_contains_required_arrays(tmp_path) -> None:
    config = ExperimentConfig(
        name="round_trip",
        experiment="free_packet",
        grid=GridConfig(n=24, length=8.0),
        wave_packet=WavePacketConfig(center=(-1.5, 0.0), sigma=0.5, momentum=(2.0, 0.0)),
        potential=PotentialConfig(kind="free"),
        solver=SolverConfig(dt=0.004, tf=0.04, frame_interval=0.02),
    )
    result = run_experiment(config)
    run_path = save_run(result, tmp_path)
    loaded = load_run(run_path)

    assert loaded["psi_frames"].ndim == 3
    assert loaded["potential"].shape == (24, 24)
    assert loaded["mask_transmission"].shape == (24, 24)
    assert loaded["metadata"]["backend"] == "numpy"

    metrics = json.loads((tmp_path / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["experiment"] == "free_packet"
