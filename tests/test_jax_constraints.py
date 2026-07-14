from __future__ import annotations

import math

import numpy as np
import pytest


jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from quantum_dynamics_lab.constraints import parameterize_phase, phase_total_variation
from quantum_dynamics_lab.jax_constraints import (
    parameterize_phase_jax,
    phase_total_variation_jax,
    quantize_phase_jax,
)


def test_jax_constraint_composition_matches_numpy_and_jits() -> None:
    rng = np.random.default_rng(20260714)
    parameters = rng.standard_normal((3, 24, 32))
    expected = parameterize_phase(
        parameters,
        bounds=(-0.75 * math.pi, 0.9 * math.pi),
        smoothing_passes=2,
        quantization_levels=16,
    )
    compiled = jax.jit(
        lambda values: parameterize_phase_jax(
            values,
            bounds=(-0.75 * math.pi, 0.9 * math.pi),
            smoothing_passes=2,
            quantization_levels=16,
            straight_through=False,
        )
    )
    actual = compiled(jnp.asarray(parameters))

    np.testing.assert_allclose(np.asarray(actual), expected, atol=2e-15, rtol=0.0)


def test_straight_through_quantization_has_hard_value_and_useful_gradient() -> None:
    parameters = jnp.linspace(-1.0, 1.0, 64).reshape(8, 8)
    hard = quantize_phase_jax(parameters, 8, straight_through=False)
    straight_through = quantize_phase_jax(parameters, 8, straight_through=True)
    gradient = jax.grad(
        lambda values: jnp.sum(
            quantize_phase_jax(values, 8, straight_through=True)
        )
    )(parameters)

    np.testing.assert_array_equal(np.asarray(straight_through), np.asarray(hard))
    np.testing.assert_allclose(np.asarray(gradient), 1.0, atol=0.0, rtol=0.0)


def test_constraint_directional_derivative_matches_centered_difference() -> None:
    rng = np.random.default_rng(4242)
    parameters = jnp.asarray(rng.standard_normal((2, 16, 16)))
    direction_values = rng.standard_normal(parameters.shape)
    direction_values /= np.linalg.norm(direction_values)
    direction = jnp.asarray(direction_values)

    def loss(values):
        phase = parameterize_phase_jax(values, smoothing_passes=2)
        return jnp.mean(jnp.sin(phase)) + 0.03 * phase_total_variation_jax(phase)

    gradient = jax.grad(loss)(parameters)
    autodiff = float(jnp.vdot(gradient, direction))
    epsilon = 1.0e-5
    finite_difference = float(
        (loss(parameters + epsilon * direction) - loss(parameters - epsilon * direction))
        / (2.0 * epsilon)
    )
    relative_error = abs(autodiff - finite_difference) / max(
        abs(autodiff),
        abs(finite_difference),
        1e-12,
    )

    assert abs(autodiff) > 1e-8
    assert relative_error < 1e-6
    phase = np.asarray(parameterize_phase_jax(parameters, smoothing_passes=2))
    assert float(phase_total_variation_jax(phase)) == pytest.approx(
        phase_total_variation(phase),
        abs=2e-15,
    )
