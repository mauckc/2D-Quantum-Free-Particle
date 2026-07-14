from __future__ import annotations

import numpy as np
import pytest


jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from quantum_dynamics_lab.jax_propagation import (
    JaxPropagationBackend,
    apply_phase_jax,
    fresnel_transfer_jax,
    multiplane_scan_jax,
    propagate_fresnel_jax,
)
from quantum_dynamics_lab.optics import gaussian_mode, make_optical_grid
from quantum_dynamics_lab.optics_propagation import (
    fresnel_transfer_function,
    propagate_fresnel,
)
from quantum_dynamics_lab.propagation import (
    prepare_split_step_operators,
    propagate_split_step,
)


@pytest.mark.parametrize(
    ("shape", "extent"),
    [((64, 64), (4.0e-3, 4.0e-3)), ((48, 64), (3.0e-3, 4.0e-3))],
)
def test_jax_fresnel_matches_numpy_x64(shape, extent) -> None:
    grid = make_optical_grid(shape, extent, 1064.0e-9)
    field = gaussian_mode(grid, 0.45e-3, center=(0.1e-3, -0.15e-3))
    distance = 0.08
    numpy_transfer = fresnel_transfer_function(grid, distance)
    numpy_field = propagate_fresnel(field, grid, distance)

    jax_transfer = fresnel_transfer_jax(
        jnp.asarray(grid.FX),
        jnp.asarray(grid.FY),
        grid.wavelength,
        distance,
    )
    jax_field = propagate_fresnel_jax(jnp.asarray(field), jax_transfer)

    assert jax.config.x64_enabled
    np.testing.assert_allclose(
        np.asarray(jax_transfer),
        numpy_transfer,
        atol=2e-10,
        rtol=2e-10,
    )
    np.testing.assert_allclose(
        np.asarray(jax_field),
        numpy_field,
        atol=2e-10,
        rtol=2e-10,
    )


def test_jax_phase_and_generic_kernel_match_numpy() -> None:
    grid = make_optical_grid(32, 3.0e-3, 633.0e-9)
    field = gaussian_mode(grid, 0.4e-3)
    phase = 0.7 * np.sin(2.0 * np.pi * grid.X / 0.8e-3)
    expected_phase = field * np.exp(1j * phase)
    backend = JaxPropagationBackend()
    operators = prepare_split_step_operators(
        -0.2j * grid.X,
        -0.1j * (grid.FX**2 + grid.FY**2),
        0.03,
        backend=backend,
    )
    result = propagate_split_step(jnp.asarray(field), operators, backend=backend)

    np.testing.assert_allclose(
        np.asarray(apply_phase_jax(field, phase)),
        expected_phase,
        atol=2e-14,
        rtol=2e-14,
    )
    assert result.dtype == jnp.complex128
    assert np.all(np.isfinite(np.asarray(result)))


def test_multiplane_scan_is_jittable_and_matches_numpy_loop() -> None:
    grid = make_optical_grid(48, 4.0e-3, 1064.0e-9)
    field = gaussian_mode(grid, 0.5e-3)
    phases = np.stack(
        [
            0.3 * np.sin(2.0 * np.pi * grid.X / 1.2e-3),
            -0.2 * np.cos(2.0 * np.pi * grid.Y / 1.0e-3),
            0.15 * np.sin(2.0 * np.pi * (grid.X + grid.Y) / 1.4e-3),
        ]
    )
    distances = (0.04, 0.05, 0.03)
    numpy_state = field
    transfers = []
    for phase, distance in zip(phases, distances, strict=True):
        numpy_state = propagate_fresnel(
            numpy_state * np.exp(1j * phase),
            grid,
            distance,
        )
        transfers.append(fresnel_transfer_function(grid, distance))

    compiled = jax.jit(multiplane_scan_jax)
    final_field, history = compiled(
        jnp.asarray(field),
        jnp.asarray(phases),
        jnp.asarray(np.stack(transfers)),
    )
    final_field.block_until_ready()

    assert jax.default_backend() == "cpu"
    assert history.shape == (3, *grid.shape)
    np.testing.assert_allclose(
        np.asarray(final_field),
        numpy_state,
        atol=3e-10,
        rtol=3e-10,
    )
