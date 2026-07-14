"""Validated NumPy reference propagators for scalar free-space optics."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from .optics import ComplexArray, FloatArray, OpticalGrid, as_complex_field

__all__ = [
    "angular_spectrum_bandlimit",
    "angular_spectrum_retained_power",
    "angular_spectrum_transfer_function",
    "fresnel_transfer_function",
    "propagate_angular_spectrum",
    "propagate_fresnel",
]


def fresnel_transfer_function(
    grid: OpticalGrid,
    distance: float,
) -> ComplexArray:
    """Return the paraxial transfer function including carrier phase."""

    distance = _finite_distance(distance)
    radial_frequency_squared = grid.FX**2 + grid.FY**2
    phase = (
        grid.wave_number * distance
        - math.pi * grid.wavelength * distance * radial_frequency_squared
    )
    return np.asarray(np.exp(1j * phase), dtype=np.complex128)


def propagate_fresnel(
    field: ArrayLike,
    grid: OpticalGrid,
    distance: float,
) -> ComplexArray:
    """Propagate a sampled field with the paraxial Fresnel transfer method."""

    values = as_complex_field(field, grid)
    transfer = fresnel_transfer_function(grid, distance)
    return np.asarray(
        np.fft.ifft2(np.fft.fft2(values) * transfer),
        dtype=np.complex128,
    )


def angular_spectrum_bandlimit(
    grid: OpticalGrid,
    distance: float,
) -> FloatArray:
    """Return the propagating, alias-safe rectangular BLAS support mask.

    The distance-dependent limits follow the finite-window band-limited angular
    spectrum condition. Nyquist limits are implicit in the sampled FFT grid.
    """

    distance = abs(_finite_distance(distance))
    if distance == 0.0:
        return np.ones(grid.shape, dtype=np.float64)
    ly, lx = grid.extent
    fx_limit = 1.0 / (
        grid.wavelength * math.sqrt(1.0 + (2.0 * distance / lx) ** 2)
    )
    fy_limit = 1.0 / (
        grid.wavelength * math.sqrt(1.0 + (2.0 * distance / ly) ** 2)
    )
    propagating = (
        grid.wavelength**2 * (grid.FX**2 + grid.FY**2)
        <= 1.0 + 8.0 * np.finfo(np.float64).eps
    )
    return np.asarray(
        propagating
        & (np.abs(grid.FX) <= fx_limit)
        & (np.abs(grid.FY) <= fy_limit),
        dtype=np.float64,
    )


def angular_spectrum_transfer_function(
    grid: OpticalGrid,
    distance: float,
    *,
    bandlimit: bool = True,
) -> ComplexArray:
    """Return the exact propagating-wave angular-spectrum transfer function."""

    distance = _finite_distance(distance)
    if distance == 0.0:
        return np.ones(grid.shape, dtype=np.complex128)
    normalized_radial_frequency = grid.wavelength**2 * (
        grid.FX**2 + grid.FY**2
    )
    propagating = normalized_radial_frequency <= 1.0
    longitudinal_factor = np.sqrt(
        np.clip(1.0 - normalized_radial_frequency, 0.0, None)
    )
    transfer = np.zeros(grid.shape, dtype=np.complex128)
    transfer[propagating] = np.exp(
        1j * grid.wave_number * distance * longitudinal_factor[propagating]
    )
    if bandlimit:
        transfer *= angular_spectrum_bandlimit(grid, distance)
    return transfer


def propagate_angular_spectrum(
    field: ArrayLike,
    grid: OpticalGrid,
    distance: float,
    *,
    bandlimit: bool = True,
) -> ComplexArray:
    """Propagate with exact scalar dispersion and optional BLAS filtering."""

    values = as_complex_field(field, grid)
    transfer = angular_spectrum_transfer_function(
        grid,
        distance,
        bandlimit=bandlimit,
    )
    return np.asarray(
        np.fft.ifft2(np.fft.fft2(values) * transfer),
        dtype=np.complex128,
    )


def angular_spectrum_retained_power(
    field: ArrayLike,
    grid: OpticalGrid,
    distance: float,
) -> float:
    """Return the fraction of sampled spectral power retained by BLAS support."""

    values = as_complex_field(field, grid)
    spectral_power = np.abs(np.fft.fft2(values)) ** 2
    total = float(np.sum(spectral_power, dtype=np.float64))
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("field spectral power must be finite and positive")
    retained = float(
        np.sum(
            spectral_power * angular_spectrum_bandlimit(grid, distance),
            dtype=np.float64,
        )
    )
    return retained / total


def _finite_distance(distance: float) -> float:
    converted = float(distance)
    if not math.isfinite(converted):
        raise ValueError("distance must be finite")
    return converted
