from __future__ import annotations

import math

import numpy as np
import pytest

from quantum_dynamics_lab.objectives import (
    intensity_similarity,
    mode_overlap,
    mode_power_efficiency,
)
from quantum_dynamics_lab.optics import (
    gaussian_mode,
    hermite_gaussian_mode,
    make_optical_grid,
)


def test_mode_objectives_have_expected_analytic_values() -> None:
    grid = make_optical_grid(128, 8.0e-3, 1064.0e-9)
    gaussian = gaussian_mode(grid, 0.6e-3)
    hg10 = hermite_gaussian_mode(grid, (1, 0), 0.6e-3)
    mixture = math.sqrt(0.8) * gaussian + math.sqrt(0.2) * hg10

    assert mode_overlap(gaussian, gaussian, grid.sample_area) == pytest.approx(
        1.0,
        abs=1e-14,
    )
    assert mode_overlap(gaussian, hg10, grid.sample_area) < 1e-28
    assert mode_overlap(mixture, gaussian, grid.sample_area) == pytest.approx(
        0.8,
        abs=1e-13,
    )
    assert mode_power_efficiency(
        mixture,
        gaussian,
        1.0,
        grid.sample_area,
    ) == pytest.approx(0.8, abs=1e-13)


def test_intensity_similarity_and_zero_power_behavior() -> None:
    field = np.array([[1.0, 2.0], [0.5, 0.25]], dtype=np.complex128)
    orthogonal_intensity = np.array([[0.0, 0.0], [0.0, 1.0]])

    assert intensity_similarity(field, field) == pytest.approx(1.0, abs=1e-15)
    assert 0.0 < intensity_similarity(field, orthogonal_intensity) < 1.0
    assert math.isnan(mode_overlap(np.zeros_like(field), field))
    assert math.isnan(intensity_similarity(np.zeros_like(field), field))
    assert math.isnan(mode_power_efficiency(field, field, 0.0))


def test_objectives_validate_shape_area_and_finite_values() -> None:
    with pytest.raises(ValueError, match="matching two-dimensional"):
        mode_overlap(np.zeros((2, 2)), np.zeros((3, 3)))
    with pytest.raises(ValueError, match="sample_area"):
        intensity_similarity(np.ones((2, 2)), np.ones((2, 2)), 0.0)
    with pytest.raises(ValueError, match="finite"):
        mode_overlap(np.full((2, 2), np.nan), np.ones((2, 2)))
