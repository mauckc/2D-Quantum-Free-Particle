"""Regenerate the small M1 scalar-optics reference metrics."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from quantum_dynamics_lab.optics import (
    apply_thin_element,
    field_power,
    gaussian_mode,
    make_optical_grid,
    thin_lens,
)
from quantum_dynamics_lab.optics_propagation import (
    angular_spectrum_retained_power,
    propagate_angular_spectrum,
    propagate_fresnel,
)


REFERENCE_PATH = Path(__file__).with_name("optics_m1.json")


def _relative_l2(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))


def _analytic_gaussian(grid, waist: float, distance: float) -> np.ndarray:
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


def _beam_waist(field: np.ndarray, grid) -> float:
    intensity = np.abs(field) ** 2
    variance = (
        float(np.sum(intensity * grid.X**2, dtype=np.float64))
        * grid.sample_area
        / field_power(field, grid)
    )
    return 2.0 * math.sqrt(variance)


def generate_metrics() -> dict[str, Any]:
    wavelength = 1064.0e-9
    waist = 0.45e-3
    distance = 0.18
    gaussian_errors = []
    for n, extent in ((64, 3.0e-3), (128, 6.0e-3)):
        grid = make_optical_grid(n, extent, wavelength)
        propagated = propagate_fresnel(gaussian_mode(grid, waist), grid, distance)
        gaussian_errors.append(
            _relative_l2(propagated, _analytic_gaussian(grid, waist, distance))
        )

    grid = make_optical_grid(192, 6.0e-3, wavelength)
    field = gaussian_mode(grid, 0.55e-3)
    fresnel = propagate_fresnel(field, grid, 0.12)
    angular = propagate_angular_spectrum(field, grid, 0.12)
    reconstructed_fresnel = propagate_fresnel(fresnel, grid, -0.12)
    reconstructed_angular = propagate_angular_spectrum(angular, grid, -0.12)

    focus_grid = make_optical_grid(256, 4.0e-3, wavelength)
    input_waist = 0.55e-3
    focal_length = 0.18
    focused = propagate_fresnel(
        apply_thin_element(
            gaussian_mode(focus_grid, input_waist),
            thin_lens(focus_grid, focal_length),
            focus_grid,
        ),
        focus_grid,
        focal_length,
    )
    predicted_focus_waist = wavelength * focal_length / (math.pi * input_waist)

    return {
        "schema_version": 1,
        "scope": "monochromatic coherent scalar free-space optics",
        "configuration": {
            "wavelength_m": wavelength,
            "gaussian_waist_m": waist,
            "gaussian_distance_m": distance,
            "cross_model_grid": 192,
            "cross_model_extent_m": 6.0e-3,
            "cross_model_waist_m": 0.55e-3,
            "cross_model_distance_m": 0.12,
        },
        "metrics": {
            "gaussian_relative_l2_coarse": gaussian_errors[0],
            "gaussian_relative_l2_fine": gaussian_errors[1],
            "gaussian_convergence_ratio": gaussian_errors[1] / gaussian_errors[0],
            "focal_waist_relative_error": abs(
                _beam_waist(focused, focus_grid) / predicted_focus_waist - 1.0
            ),
            "fresnel_power_relative_drift": abs(field_power(fresnel, grid) - 1.0),
            "fresnel_reverse_max_error": float(
                np.max(np.abs(reconstructed_fresnel - field))
            ),
            "angular_retained_spectral_power": angular_spectrum_retained_power(
                field, grid, 0.12
            ),
            "angular_power_relative_drift": abs(field_power(angular, grid) - 1.0),
            "angular_reverse_max_error": float(
                np.max(np.abs(reconstructed_angular - field))
            ),
            "cross_model_complex_relative_l2": _relative_l2(angular, fresnel),
            "cross_model_intensity_relative_l2": _relative_l2(
                np.abs(angular) ** 2,
                np.abs(fresnel) ** 2,
            ),
        },
        "declared_limits": {
            "gaussian_relative_l2_fine_max": 1.0e-9,
            "gaussian_convergence_ratio_max": 1.0e-4,
            "focal_waist_relative_error_max": 1.0e-2,
            "power_relative_drift_max": 1.0e-13,
            "reverse_max_error_max": 1.0e-10,
            "cross_model_relative_l2_max": 5.0e-6,
            "minimum_retained_spectral_power": 1.0 - 1.0e-13,
        },
    }


def main() -> None:
    REFERENCE_PATH.write_text(
        json.dumps(generate_metrics(), indent=2) + "\n",
        encoding="utf-8",
    )
    print(REFERENCE_PATH)


if __name__ == "__main__":
    main()
