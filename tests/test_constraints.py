from __future__ import annotations

import math

import numpy as np
import pytest

from quantum_dynamics_lab.constraints import (
    bounded_phase,
    parameterize_phase,
    phase_total_variation,
    quantize_phase,
    smooth_periodic,
)


def test_bounded_phase_and_composition_stay_inside_limits() -> None:
    parameters = np.linspace(-20.0, 20.0, 256).reshape(16, 16)
    bounds = (-0.25 * math.pi, 0.75 * math.pi)

    phase = bounded_phase(parameters, bounds)
    composed = parameterize_phase(parameters, bounds=bounds, smoothing_passes=2)

    for values in (phase, composed):
        assert np.min(values) >= bounds[0]
        assert np.max(values) <= bounds[1]
        assert np.all(np.isfinite(values))


def test_smoothing_reduces_deterministic_total_variation() -> None:
    checkerboard = (-1.0) ** sum(np.indices((32, 32)))
    phase = bounded_phase(checkerboard)
    smoothed = bounded_phase(smooth_periodic(checkerboard, passes=3))

    assert phase_total_variation(smoothed) < 0.1 * phase_total_variation(phase)


def test_wrapped_tv_is_invariant_to_full_phase_cycles() -> None:
    rng = np.random.default_rng(42)
    phase = rng.uniform(-math.pi, math.pi, (16, 16))
    cycles = 2.0 * math.pi * rng.integers(-3, 4, phase.shape)

    assert phase_total_variation(phase + cycles) == pytest.approx(
        phase_total_variation(phase),
        abs=2e-15,
    )


def test_quantization_has_exact_requested_periodic_levels() -> None:
    phase = np.linspace(-math.pi, math.pi, 4096, endpoint=False).reshape(64, 64)
    quantized = quantize_phase(phase, 8)
    expected = -math.pi + np.arange(8) * (2.0 * math.pi / 8)

    np.testing.assert_allclose(np.unique(quantized), expected, atol=1e-15)
    np.testing.assert_allclose(
        parameterize_phase(phase, quantization_levels=8),
        quantize_phase(bounded_phase(phase), 8),
        atol=0.0,
    )


@pytest.mark.parametrize(
    "call",
    [
        lambda: bounded_phase(np.ones((2, 2)), (1.0, 1.0)),
        lambda: smooth_periodic(np.ones((2, 2)), -1),
        lambda: quantize_phase(np.ones((2, 2)), 1),
    ],
)
def test_invalid_constraint_configuration_fails(call) -> None:
    with pytest.raises(ValueError):
        call()
