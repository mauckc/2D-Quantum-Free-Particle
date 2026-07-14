from __future__ import annotations

import numpy as np
import pytest

from quantum_dynamics_lab.optics import (
    apply_thin_element,
    as_complex_field,
    circular_aperture,
    field_power,
    gaussian_mode,
    hermite_gaussian_mode,
    make_optical_grid,
    phase_only_mask,
    rectangular_aperture,
    thin_lens,
)


def test_optical_grid_uses_physical_sampling_and_read_only_coordinates() -> None:
    grid = make_optical_grid((48, 64), (3.0e-3, 4.0e-3), 633.0e-9)

    assert grid.shape == (48, 64)
    assert grid.dy == pytest.approx(3.0e-3 / 48)
    assert grid.dx == pytest.approx(4.0e-3 / 64)
    assert grid.wave_number == pytest.approx(2.0 * np.pi / 633.0e-9)
    assert grid.X.shape == grid.shape
    assert grid.FX.shape == grid.shape
    assert not grid.X.flags.writeable


@pytest.mark.parametrize(
    ("shape", "extent", "wavelength"),
    [
        (1, 1.0, 1.0),
        (16, 0.0, 1.0),
        (16, 1.0, -1.0),
        ((16, 0), (1.0, 1.0), 1.0),
    ],
)
def test_optical_grid_rejects_nonphysical_inputs(shape, extent, wavelength) -> None:
    with pytest.raises(ValueError):
        make_optical_grid(shape, extent, wavelength)


def test_arbitrary_complex_field_is_copied_without_normalization() -> None:
    grid = make_optical_grid(8, 2.0e-3, 532.0e-9)
    values = np.full(grid.shape, 2.0 + 3.0j)

    field = as_complex_field(values, grid)
    field[0, 0] = 9.0

    assert values[0, 0] == 2.0 + 3.0j
    assert field_power(as_complex_field(values, grid), grid) == pytest.approx(
        13.0 * 4.0e-6
    )


def test_gaussian_and_hg_modes_are_normalized_and_orthogonal() -> None:
    grid = make_optical_grid(128, 8.0e-3, 1064.0e-9)
    gaussian = gaussian_mode(grid, 0.7e-3)
    hg10 = hermite_gaussian_mode(grid, (1, 0), 0.7e-3)
    overlap = np.vdot(gaussian, hg10) * grid.sample_area

    assert field_power(gaussian, grid) == pytest.approx(1.0, abs=1e-13)
    assert field_power(hg10, grid) == pytest.approx(1.0, abs=1e-13)
    assert abs(overlap) < 1e-14


def test_apertures_are_binary_and_cannot_increase_power() -> None:
    grid = make_optical_grid(64, 4.0e-3, 532.0e-9)
    field = gaussian_mode(grid, 0.8e-3)
    circular = circular_aperture(grid, 0.6e-3)
    rectangular = rectangular_aperture(grid, (1.0e-3, 1.5e-3))

    for aperture in (circular, rectangular):
        assert set(np.unique(aperture)) <= {0.0, 1.0}
        truncated = apply_thin_element(field, aperture, grid)
        assert field_power(truncated, grid) <= field_power(field, grid)


def test_lens_matches_paraxial_phase_and_preserves_power() -> None:
    grid = make_optical_grid(64, 3.0e-3, 633.0e-9)
    field = gaussian_mode(grid, 0.5e-3)
    focal_length = 0.15
    lens = thin_lens(grid, focal_length)
    expected = np.exp(
        -0.5j
        * grid.wave_number
        * (grid.X**2 + grid.Y**2)
        / focal_length
    )
    output = apply_thin_element(field, lens, grid)

    np.testing.assert_allclose(lens, expected, atol=1e-14, rtol=1e-14)
    assert field_power(output, grid) == pytest.approx(1.0, abs=1e-13)


def test_phase_mask_is_bounded_unit_magnitude_and_pure() -> None:
    grid = make_optical_grid(32, 2.0e-3, 633.0e-9)
    phase = np.linspace(-np.pi, np.pi, np.prod(grid.shape)).reshape(grid.shape)
    before = phase.copy()
    mask = phase_only_mask(phase, grid)

    np.testing.assert_array_equal(phase, before)
    np.testing.assert_allclose(np.abs(mask), 1.0, atol=1e-15)
    with pytest.raises(ValueError, match="within bounds"):
        phase_only_mask(phase + 0.1, grid)


def test_thin_element_rejects_gain_and_shape_mismatch() -> None:
    grid = make_optical_grid(8, 1.0e-3, 633.0e-9)
    with pytest.raises(ValueError, match="cannot exceed one"):
        apply_thin_element(np.ones(grid.shape), np.full(grid.shape, 1.01), grid)
    with pytest.raises(ValueError, match="shape"):
        as_complex_field(np.ones((4, 4)), grid)
