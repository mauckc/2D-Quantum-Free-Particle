from __future__ import annotations

import math

import numpy as np
import pytest

from quantum_dynamics_lab.optics import (
    apply_thin_element,
    field_power,
    gaussian_mode,
    make_optical_grid,
    thin_lens,
)
from quantum_dynamics_lab.optics_propagation import (
    angular_spectrum_bandlimit,
    angular_spectrum_retained_power,
    propagate_angular_spectrum,
    propagate_fresnel,
)


def _analytic_gaussian(grid, waist: float, distance: float):
    envelope = np.exp(-(grid.X**2 + grid.Y**2) / waist**2)
    amplitude = 1.0 / math.sqrt(field_power(envelope, grid))
    rayleigh_range = math.pi * waist**2 / grid.wavelength
    q = 1.0 + 1j * distance / rayleigh_range
    return np.asarray(
        amplitude
        * np.exp(1j * grid.wave_number * distance)
        * np.exp(-(grid.X**2 + grid.Y**2) / (waist**2 * q))
        / q,
        dtype=np.complex128,
    )


def _relative_l2(actual, expected) -> float:
    return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))


def _beam_waist(field, grid) -> float:
    intensity = np.abs(field) ** 2
    power = field_power(field, grid)
    x_variance = (
        float(np.sum(intensity * grid.X**2, dtype=np.float64))
        * grid.sample_area
        / power
    )
    return 2.0 * math.sqrt(x_variance)


def test_fresnel_matches_analytic_gaussian_and_converges() -> None:
    wavelength = 1064.0e-9
    waist = 0.45e-3
    distance = 0.18
    errors = []

    for n, extent in ((64, 3.0e-3), (128, 6.0e-3)):
        grid = make_optical_grid(n, extent, wavelength)
        propagated = propagate_fresnel(
            gaussian_mode(grid, waist),
            grid,
            distance,
        )
        errors.append(
            _relative_l2(propagated, _analytic_gaussian(grid, waist, distance))
        )

    assert errors[1] < errors[0] * 1e-4
    assert errors[1] < 1e-9


def test_fresnel_is_power_conserving_and_reversible() -> None:
    grid = make_optical_grid((48, 64), (4.0e-3, 5.0e-3), 633.0e-9)
    rng = np.random.default_rng(34)
    field = rng.standard_normal(grid.shape) + 1j * rng.standard_normal(grid.shape)
    before = field_power(field, grid)

    propagated = propagate_fresnel(field, grid, 0.07)
    reconstructed = propagate_fresnel(propagated, grid, -0.07)

    assert abs(field_power(propagated, grid) - before) / before < 1e-14
    assert np.max(np.abs(reconstructed - field)) < 3e-15


def test_lens_focus_matches_predicted_gaussian_waist() -> None:
    wavelength = 1064.0e-9
    input_waist = 0.55e-3
    focal_length = 0.18
    grid = make_optical_grid(256, 4.0e-3, wavelength)
    field = gaussian_mode(grid, input_waist)
    focused = propagate_fresnel(
        apply_thin_element(field, thin_lens(grid, focal_length), grid),
        grid,
        focal_length,
    )
    expected_waist = wavelength * focal_length / (math.pi * input_waist)

    assert abs(_beam_waist(focused, grid) / expected_waist - 1.0) < 0.01


def test_angular_spectrum_is_conservative_and_reversible_in_retained_band() -> None:
    grid = make_optical_grid(128, 5.0e-3, 1064.0e-9)
    field = gaussian_mode(grid, 0.5e-3)
    retained = angular_spectrum_retained_power(field, grid, 0.12)

    propagated = propagate_angular_spectrum(field, grid, 0.12)
    reconstructed = propagate_angular_spectrum(propagated, grid, -0.12)

    assert 1.0 - retained < 1e-14
    assert abs(field_power(propagated, grid) - 1.0) < 1e-13
    assert np.max(np.abs(reconstructed - field)) < 1e-11


def test_bandlimit_rejects_alias_prone_spectrum() -> None:
    grid = make_optical_grid(64, 1.0e-3, 532.0e-9)
    mask = angular_spectrum_bandlimit(grid, 2.0)
    checkerboard = (-1.0) ** (
        np.indices(grid.shape)[0] + np.indices(grid.shape)[1]
    )

    assert np.count_nonzero(mask) < mask.size
    assert angular_spectrum_retained_power(checkerboard, grid, 2.0) < 1e-15
    rejected = propagate_angular_spectrum(checkerboard, grid, 2.0)
    assert field_power(rejected, grid) < 1e-20


def test_models_agree_in_declared_paraxial_regime() -> None:
    grid = make_optical_grid(192, 6.0e-3, 1064.0e-9)
    field = gaussian_mode(grid, 0.55e-3)
    distance = 0.12

    fresnel = propagate_fresnel(field, grid, distance)
    angular = propagate_angular_spectrum(field, grid, distance)

    assert angular_spectrum_retained_power(field, grid, distance) > 1.0 - 1e-14
    assert _relative_l2(angular, fresnel) < 5e-6
    intensity_error = _relative_l2(np.abs(angular) ** 2, np.abs(fresnel) ** 2)
    assert intensity_error < 5e-6


def test_propagators_reject_nonfinite_distance_and_preserve_inputs() -> None:
    grid = make_optical_grid(16, 1.0e-3, 633.0e-9)
    field = gaussian_mode(grid, 0.2e-3)
    before = field.copy()
    with pytest.raises(ValueError, match="distance"):
        propagate_fresnel(field, grid, float("nan"))
    with pytest.raises(ValueError, match="distance"):
        propagate_angular_spectrum(field, grid, float("inf"))
    np.testing.assert_array_equal(field, before)
