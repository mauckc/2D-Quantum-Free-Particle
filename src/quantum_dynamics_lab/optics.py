"""Pure NumPy foundations for coherent monochromatic scalar optics."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]

__all__ = [
    "OpticalGrid",
    "apply_thin_element",
    "as_complex_field",
    "circular_aperture",
    "field_power",
    "gaussian_mode",
    "hermite_gaussian_mode",
    "make_optical_grid",
    "normalize_field",
    "phase_only_mask",
    "rectangular_aperture",
    "thin_lens",
]


@dataclass(frozen=True)
class OpticalGrid:
    """Sampled transverse plane for a scalar field.

    ``shape`` and ``extent`` are ordered as ``(y, x)``. Spatial-frequency
    arrays use NumPy FFT ordering so they can multiply ``fft2`` output directly.
    """

    shape: tuple[int, int]
    extent: tuple[float, float]
    wavelength: float
    dy: float
    dx: float
    y: FloatArray
    x: FloatArray
    Y: FloatArray
    X: FloatArray
    fy: FloatArray
    fx: FloatArray
    FY: FloatArray
    FX: FloatArray

    @property
    def wave_number(self) -> float:
        return 2.0 * math.pi / self.wavelength

    @property
    def sample_area(self) -> float:
        return self.dy * self.dx


def make_optical_grid(
    shape: int | tuple[int, int],
    extent: float | tuple[float, float],
    wavelength: float,
) -> OpticalGrid:
    """Construct a centered periodic grid using consistent physical units."""

    ny, nx = _pair_of_ints(shape, "shape")
    ly, lx = _pair_of_positive_floats(extent, "extent")
    if not math.isfinite(wavelength) or wavelength <= 0.0:
        raise ValueError("wavelength must be finite and positive")

    dy = ly / ny
    dx = lx / nx
    y = (np.arange(ny, dtype=np.float64) - ny // 2) * dy
    x = (np.arange(nx, dtype=np.float64) - nx // 2) * dx
    fy = np.fft.fftfreq(ny, d=dy).astype(np.float64)
    fx = np.fft.fftfreq(nx, d=dx).astype(np.float64)
    X, Y = np.meshgrid(x, y, indexing="xy")
    FX, FY = np.meshgrid(fx, fy, indexing="xy")
    arrays = tuple(_read_only(values) for values in (y, x, Y, X, fy, fx, FY, FX))
    return OpticalGrid(
        shape=(ny, nx),
        extent=(ly, lx),
        wavelength=float(wavelength),
        dy=dy,
        dx=dx,
        y=arrays[0],
        x=arrays[1],
        Y=arrays[2],
        X=arrays[3],
        fy=arrays[4],
        fx=arrays[5],
        FY=arrays[6],
        FX=arrays[7],
    )


def as_complex_field(values: ArrayLike, grid: OpticalGrid) -> ComplexArray:
    """Copy an arbitrary complex field without changing its amplitude."""

    field = np.array(values, dtype=np.complex128, copy=True)
    _validate_shape(field, grid, "field")
    if not np.all(np.isfinite(field)):
        raise ValueError("field must contain only finite values")
    return field


def field_power(field: ArrayLike, grid: OpticalGrid) -> float:
    """Integrate scalar-field power over the sampled transverse plane."""

    values = np.asarray(field)
    _validate_shape(values, grid, "field")
    return float(np.sum(np.abs(values) ** 2, dtype=np.float64) * grid.sample_area)


def normalize_field(field: ArrayLike, grid: OpticalGrid) -> ComplexArray:
    """Return a unit-power complex128 copy of ``field``."""

    values = as_complex_field(field, grid)
    power = field_power(values, grid)
    if not math.isfinite(power) or power <= 0.0:
        raise ValueError("field power must be finite and positive")
    return np.asarray(values / math.sqrt(power), dtype=np.complex128)


def gaussian_mode(
    grid: OpticalGrid,
    waist: float,
    *,
    center: tuple[float, float] = (0.0, 0.0),
) -> ComplexArray:
    """Return a unit-power waist-plane Gaussian field.

    ``waist`` is the field-amplitude 1/e radius and ``center`` is ``(y, x)``.
    """

    waist = _positive_finite(waist, "waist")
    cy, cx = _finite_pair(center, "center")
    radius_squared = (grid.X - cx) ** 2 + (grid.Y - cy) ** 2
    return normalize_field(np.exp(-radius_squared / waist**2), grid)


def hermite_gaussian_mode(
    grid: OpticalGrid,
    order: tuple[int, int],
    waist: float,
    *,
    center: tuple[float, float] = (0.0, 0.0),
) -> ComplexArray:
    """Return a unit-power waist-plane Hermite-Gaussian ``(x, y)`` mode."""

    mx, my = _pair_of_nonnegative_ints(order, "order")
    waist = _positive_finite(waist, "waist")
    cy, cx = _finite_pair(center, "center")
    scaled_y = math.sqrt(2.0) * (grid.Y - cy) / waist
    scaled_x = math.sqrt(2.0) * (grid.X - cx) / waist
    coefficients_y = np.zeros(my + 1, dtype=np.float64)
    coefficients_x = np.zeros(mx + 1, dtype=np.float64)
    coefficients_y[my] = 1.0
    coefficients_x[mx] = 1.0
    envelope = np.exp(
        -((grid.X - cx) ** 2 + (grid.Y - cy) ** 2) / waist**2
    )
    field = (
        np.polynomial.hermite.hermval(scaled_y, coefficients_y)
        * np.polynomial.hermite.hermval(scaled_x, coefficients_x)
        * envelope
    )
    return normalize_field(field, grid)


def circular_aperture(
    grid: OpticalGrid,
    radius: float,
    *,
    center: tuple[float, float] = (0.0, 0.0),
) -> FloatArray:
    """Return the binary transmission of a circular aperture."""

    radius = _positive_finite(radius, "radius")
    cy, cx = _finite_pair(center, "center")
    return np.asarray(
        (grid.X - cx) ** 2 + (grid.Y - cy) ** 2 <= radius**2,
        dtype=np.float64,
    )


def rectangular_aperture(
    grid: OpticalGrid,
    size: float | tuple[float, float],
    *,
    center: tuple[float, float] = (0.0, 0.0),
) -> FloatArray:
    """Return the binary transmission of a centered ``(height, width)`` aperture."""

    height, width = _pair_of_positive_floats(size, "size")
    cy, cx = _finite_pair(center, "center")
    return np.asarray(
        (np.abs(grid.Y - cy) <= height / 2.0)
        & (np.abs(grid.X - cx) <= width / 2.0),
        dtype=np.float64,
    )


def thin_lens(
    grid: OpticalGrid,
    focal_length: float,
    *,
    center: tuple[float, float] = (0.0, 0.0),
) -> ComplexArray:
    """Return paraxial thin-lens transmission for the exp(+i k z) convention."""

    if not math.isfinite(focal_length) or focal_length == 0.0:
        raise ValueError("focal_length must be finite and nonzero")
    cy, cx = _finite_pair(center, "center")
    phase = -grid.wave_number * (
        (grid.X - cx) ** 2 + (grid.Y - cy) ** 2
    ) / (2.0 * focal_length)
    return np.asarray(np.exp(1j * phase), dtype=np.complex128)


def phase_only_mask(
    phase: ArrayLike,
    grid: OpticalGrid,
    *,
    bounds: tuple[float, float] = (-math.pi, math.pi),
) -> ComplexArray:
    """Convert a finite, bounded phase array to unit-magnitude transmission."""

    values = np.asarray(phase, dtype=np.float64)
    _validate_shape(values, grid, "phase")
    lower, upper = _finite_pair(bounds, "bounds")
    if lower >= upper:
        raise ValueError("phase bounds must be strictly increasing")
    if not np.all(np.isfinite(values)):
        raise ValueError("phase must contain only finite values")
    tolerance = 16.0 * np.finfo(np.float64).eps
    if np.any(values < lower - tolerance) or np.any(values > upper + tolerance):
        raise ValueError("phase values must lie within bounds")
    return np.asarray(np.exp(1j * values), dtype=np.complex128)


def apply_thin_element(
    field: ArrayLike,
    transmission: ArrayLike,
    grid: OpticalGrid,
) -> ComplexArray:
    """Apply a sampled thin element without mutating either input."""

    values = as_complex_field(field, grid)
    element = np.asarray(transmission, dtype=np.complex128)
    _validate_shape(element, grid, "transmission")
    if not np.all(np.isfinite(element)):
        raise ValueError("transmission must contain only finite values")
    if np.any(np.abs(element) > 1.0 + 16.0 * np.finfo(np.float64).eps):
        raise ValueError("passive transmission magnitude cannot exceed one")
    return np.asarray(values * element, dtype=np.complex128)


def _pair_of_ints(values: int | tuple[int, int], name: str) -> tuple[int, int]:
    pair = (values, values) if isinstance(values, int) else values
    if len(pair) != 2 or any(
        isinstance(value, bool) or not isinstance(value, int) for value in pair
    ):
        raise ValueError(f"{name} must contain two integers")
    if any(value < 2 for value in pair):
        raise ValueError(f"{name} entries must be at least two")
    return pair


def _pair_of_nonnegative_ints(
    values: tuple[int, int], name: str
) -> tuple[int, int]:
    if len(values) != 2 or any(
        isinstance(value, bool) or not isinstance(value, int) for value in values
    ):
        raise ValueError(f"{name} must contain two integers")
    if any(value < 0 for value in values):
        raise ValueError(f"{name} entries must be nonnegative")
    return values


def _pair_of_positive_floats(
    values: float | tuple[float, float], name: str
) -> tuple[float, float]:
    pair = (values, values) if isinstance(values, (int, float)) else values
    if len(pair) != 2:
        raise ValueError(f"{name} must contain two values")
    converted = (float(pair[0]), float(pair[1]))
    if any(not math.isfinite(value) or value <= 0.0 for value in converted):
        raise ValueError(f"{name} entries must be finite and positive")
    return converted


def _finite_pair(values: tuple[float, float], name: str) -> tuple[float, float]:
    if len(values) != 2:
        raise ValueError(f"{name} must contain two values")
    converted = (float(values[0]), float(values[1]))
    if any(not math.isfinite(value) for value in converted):
        raise ValueError(f"{name} entries must be finite")
    return converted


def _positive_finite(value: float, name: str) -> float:
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return converted


def _validate_shape(values: NDArray[np.generic], grid: OpticalGrid, name: str) -> None:
    if values.ndim != 2 or values.shape != grid.shape:
        raise ValueError(f"{name} must have shape {grid.shape}")


def _read_only(values: FloatArray) -> FloatArray:
    values.setflags(write=False)
    return values
