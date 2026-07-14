from __future__ import annotations

import numpy as np
import pytest

from quantum_dynamics_lab.backends import NumPyBackend
from quantum_dynamics_lab.config import GridConfig, SolverConfig, WavePacketConfig
from quantum_dynamics_lab.grid import integrate_probability, make_grid, norm, normalize
from quantum_dynamics_lab.propagation import (
    NumPyPropagationBackend,
    prepare_split_step_operators,
    propagate_split_step,
)
from quantum_dynamics_lab.solver import SplitStepSolver
from quantum_dynamics_lab.wavepackets import gaussian_packet


def test_plane_wave_accumulates_expected_free_phase() -> None:
    grid = make_grid(GridConfig(n=32, length=8.0))
    mode = 2
    wave_number = 2.0 * np.pi * mode / grid.length
    field = normalize(
        np.exp(1j * wave_number * grid.X).astype(np.complex128), grid
    )
    step_size = 0.01
    backend = NumPyPropagationBackend()
    operators = prepare_split_step_operators(
        np.zeros_like(grid.X, dtype=np.complex128),
        -0.5j * grid.k2,
        step_size,
        backend=backend,
    )

    evolved = propagate_split_step(field, operators, backend=backend)
    expected = field * np.exp(-0.5j * wave_number**2 * step_size)

    assert operators.position_half_step.dtype == np.complex128
    assert operators.spectral_step.dtype == np.complex128
    assert evolved.dtype == np.complex128
    assert np.max(np.abs(evolved - expected)) < 1e-12
    assert abs(norm(evolved, grid) - 1.0) < 1e-12


def test_real_potential_step_is_conservative_and_reversible() -> None:
    grid = make_grid(GridConfig(n=32, length=8.0))
    rng = np.random.default_rng(20260713)
    field = normalize(
        (
            rng.standard_normal((grid.n, grid.n))
            + 1j * rng.standard_normal((grid.n, grid.n))
        ).astype(np.complex128),
        grid,
    )
    potential = 0.05 * (grid.X**2 + 0.5 * grid.Y**2)
    position_generator = -1j * potential
    spectral_generator = -0.5j * grid.k2
    backend = NumPyPropagationBackend()
    forward = prepare_split_step_operators(
        position_generator,
        spectral_generator,
        0.004,
        backend=backend,
    )
    reverse = prepare_split_step_operators(
        position_generator,
        spectral_generator,
        -0.004,
        backend=backend,
    )

    evolved = propagate_split_step(field, forward, backend=backend)
    reconstructed = propagate_split_step(evolved, reverse, backend=backend)

    assert abs(norm(evolved, grid) - norm(field, grid)) < 1e-12
    assert np.max(np.abs(reconstructed - field)) < 1e-12


def test_free_gaussian_variance_converges_with_grid_refinement() -> None:
    sigma = 0.6
    propagation_time = 0.4
    expected_variance = sigma**2 + propagation_time**2 / (4.0 * sigma**2)
    errors = []
    backend = NumPyPropagationBackend()

    for n in (16, 24):
        grid = make_grid(GridConfig(n=n, length=12.0))
        field = gaussian_packet(
            WavePacketConfig(
                center=(0.0, 0.0),
                sigma=sigma,
                momentum=(0.0, 0.0),
            ),
            grid,
        )
        operators = prepare_split_step_operators(
            np.zeros_like(grid.X),
            -0.5j * grid.k2,
            propagation_time,
            backend=backend,
        )
        evolved = propagate_split_step(field, operators, backend=backend)
        probability = np.abs(evolved) ** 2
        mean_x = integrate_probability(probability * grid.X, grid)
        variance_x = integrate_probability(
            probability * (grid.X - mean_x) ** 2,
            grid,
        )
        errors.append(abs(variance_x - expected_variance))

    assert errors[1] < errors[0] * 1e-4
    assert errors[1] < 1e-9


@pytest.mark.parametrize("potential_kind", ["free", "static"])
def test_pure_kernel_matches_pre_refactor_solver(potential_kind: str) -> None:
    grid = make_grid(GridConfig(n=32, length=8.0))
    rng = np.random.default_rng(17)
    field = normalize(
        (
            rng.standard_normal((grid.n, grid.n))
            + 1j * rng.standard_normal((grid.n, grid.n))
        ).astype(np.complex128),
        grid,
    )
    if potential_kind == "free":
        potential = np.zeros((grid.n, grid.n), dtype=np.float64)
    else:
        potential = 0.2 * np.cos(2.0 * np.pi * grid.X / grid.length) + 0.1 * np.sin(
            2.0 * np.pi * grid.Y / grid.length
        )
    config = SolverConfig(
        dt=0.007,
        tf=0.007,
        frame_interval=0.007,
        hbar=1.3,
        mass=0.8,
        backend="numpy",
    )
    legacy = SplitStepSolver(
        grid,
        potential,
        config,
        backend=NumPyBackend(),
    ).step(field)
    backend = NumPyPropagationBackend()
    operators = prepare_split_step_operators(
        -1j * potential / config.hbar,
        -1j * config.hbar * grid.k2 / (2.0 * config.mass),
        config.dt,
        backend=backend,
    )

    propagated = propagate_split_step(field, operators, backend=backend)

    assert np.max(np.abs(propagated - legacy)) < 1e-12


def test_inputs_must_be_matching_two_dimensional_arrays() -> None:
    backend = NumPyPropagationBackend()
    with pytest.raises(ValueError, match="generators must be two-dimensional"):
        prepare_split_step_operators(
            np.zeros(4),
            np.zeros(4),
            0.1,
            backend=backend,
        )
    operators = prepare_split_step_operators(
        np.zeros((4, 4)),
        np.zeros((4, 4)),
        0.1,
        backend=backend,
    )
    with pytest.raises(ValueError, match="must have matching shapes"):
        propagate_split_step(np.zeros((8, 8)), operators, backend=backend)


def test_propagation_does_not_mutate_inputs() -> None:
    backend = NumPyPropagationBackend()
    field = np.arange(16, dtype=np.float64).reshape(4, 4).astype(np.complex128)
    operators = prepare_split_step_operators(
        -0.1j * np.ones((4, 4)),
        -0.2j * np.ones((4, 4)),
        0.05,
        backend=backend,
    )
    field_before = field.copy()
    position_before = operators.position_half_step.copy()
    spectral_before = operators.spectral_step.copy()

    result = propagate_split_step(field, operators, backend=backend)

    np.testing.assert_array_equal(field, field_before)
    np.testing.assert_array_equal(operators.position_half_step, position_before)
    np.testing.assert_array_equal(operators.spectral_step, spectral_before)
    assert not np.shares_memory(result, field)
