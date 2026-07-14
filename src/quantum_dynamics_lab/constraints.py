"""NumPy reference phase parameterization and fabrication constraints."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]

__all__ = [
    "bounded_phase",
    "parameterize_phase",
    "phase_total_variation",
    "quantize_phase",
    "smooth_periodic",
]


def bounded_phase(
    parameters: ArrayLike,
    bounds: tuple[float, float] = (-math.pi, math.pi),
) -> FloatArray:
    """Map unconstrained real variables smoothly inside phase bounds."""

    lower, upper = _validate_bounds(bounds)
    values = _finite_array(parameters)
    midpoint = 0.5 * (lower + upper)
    half_range = 0.5 * (upper - lower)
    return np.asarray(midpoint + half_range * np.tanh(values), dtype=np.float64)


def smooth_periodic(values: ArrayLike, passes: int = 1) -> FloatArray:
    """Apply repeated differentiable five-point periodic smoothing."""

    passes = _validate_passes(passes)
    smoothed = _finite_array(values).copy()
    for _ in range(passes):
        smoothed = (
            0.5 * smoothed
            + 0.125 * np.roll(smoothed, 1, axis=-2)
            + 0.125 * np.roll(smoothed, -1, axis=-2)
            + 0.125 * np.roll(smoothed, 1, axis=-1)
            + 0.125 * np.roll(smoothed, -1, axis=-1)
        )
    return np.asarray(smoothed, dtype=np.float64)


def phase_total_variation(phase: ArrayLike, epsilon: float = 1.0e-8) -> float:
    """Return mean isotropic TV using wrapped periodic phase differences."""

    values = _finite_array(phase)
    epsilon = float(epsilon)
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    dx = _wrapped_difference(np.roll(values, -1, axis=-1) - values)
    dy = _wrapped_difference(np.roll(values, -1, axis=-2) - values)
    return float(np.mean(np.sqrt(dx**2 + dy**2 + epsilon**2) - epsilon))


def quantize_phase(
    phase: ArrayLike,
    levels: int,
    bounds: tuple[float, float] = (-math.pi, math.pi),
) -> FloatArray:
    """Hard-quantize periodic phase into equally spaced levels in [min, max)."""

    levels = _validate_levels(levels)
    lower, upper = _validate_bounds(bounds)
    values = _finite_array(phase)
    span = upper - lower
    step = span / levels
    wrapped = np.mod(values - lower, span)
    indices = np.mod(np.floor(wrapped / step + 0.5), levels)
    return np.asarray(lower + indices * step, dtype=np.float64)


def parameterize_phase(
    parameters: ArrayLike,
    *,
    bounds: tuple[float, float] = (-math.pi, math.pi),
    smoothing_passes: int = 0,
    quantization_levels: int | None = None,
) -> FloatArray:
    """Compose smoothing, bounded mapping, and optional hard quantization."""

    values = smooth_periodic(parameters, smoothing_passes)
    phase = bounded_phase(values, bounds)
    if quantization_levels is not None:
        phase = quantize_phase(phase, quantization_levels, bounds)
    return phase


def _wrapped_difference(values: FloatArray) -> FloatArray:
    return np.arctan2(np.sin(values), np.cos(values))


def _finite_array(values: ArrayLike) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim < 2:
        raise ValueError("phase arrays must have at least two dimensions")
    if not np.all(np.isfinite(array)):
        raise ValueError("phase arrays must contain only finite values")
    return array


def _validate_bounds(bounds: tuple[float, float]) -> tuple[float, float]:
    if len(bounds) != 2:
        raise ValueError("bounds must contain two values")
    lower, upper = float(bounds[0]), float(bounds[1])
    if not math.isfinite(lower) or not math.isfinite(upper) or lower >= upper:
        raise ValueError("bounds must be finite and strictly increasing")
    return lower, upper


def _validate_passes(passes: int) -> int:
    if isinstance(passes, bool) or not isinstance(passes, int) or passes < 0:
        raise ValueError("smoothing passes must be a nonnegative integer")
    return passes


def _validate_levels(levels: int) -> int:
    if isinstance(levels, bool) or not isinstance(levels, int) or levels < 2:
        raise ValueError("quantization levels must be an integer of at least two")
    return levels
