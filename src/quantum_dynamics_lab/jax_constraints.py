"""Optional JAX phase parameterization and fabrication constraints."""

from __future__ import annotations

import math
from typing import Any

import jax
import jax.numpy as jnp


def bounded_phase_jax(
    parameters: Any,
    bounds: tuple[float, float] = (-math.pi, math.pi),
):
    lower, upper = _validate_bounds(bounds)
    values = jnp.asarray(parameters)
    midpoint = 0.5 * (lower + upper)
    half_range = 0.5 * (upper - lower)
    return midpoint + half_range * jnp.tanh(values)


def smooth_periodic_jax(values: Any, passes: int = 1):
    passes = _validate_passes(passes)
    smoothed = jnp.asarray(values)
    for _ in range(passes):
        smoothed = (
            0.5 * smoothed
            + 0.125 * jnp.roll(smoothed, 1, axis=-2)
            + 0.125 * jnp.roll(smoothed, -1, axis=-2)
            + 0.125 * jnp.roll(smoothed, 1, axis=-1)
            + 0.125 * jnp.roll(smoothed, -1, axis=-1)
        )
    return smoothed


def phase_total_variation_jax(phase: Any, epsilon: float = 1.0e-8):
    values = jnp.asarray(phase)
    dx = _wrapped_difference_jax(jnp.roll(values, -1, axis=-1) - values)
    dy = _wrapped_difference_jax(jnp.roll(values, -1, axis=-2) - values)
    return jnp.mean(jnp.sqrt(dx**2 + dy**2 + epsilon**2) - epsilon)


def quantize_phase_jax(
    phase: Any,
    levels: int,
    bounds: tuple[float, float] = (-math.pi, math.pi),
    *,
    straight_through: bool = False,
):
    levels = _validate_levels(levels)
    lower, upper = _validate_bounds(bounds)
    values = jnp.asarray(phase)
    span = upper - lower
    step = span / levels
    wrapped = jnp.mod(values - lower, span)
    indices = jnp.mod(jnp.floor(wrapped / step + 0.5), levels)
    quantized = lower + indices * step
    if straight_through:
        return values + jax.lax.stop_gradient(quantized - values)
    return quantized


def parameterize_phase_jax(
    parameters: Any,
    *,
    bounds: tuple[float, float] = (-math.pi, math.pi),
    smoothing_passes: int = 0,
    quantization_levels: int | None = None,
    straight_through: bool = True,
):
    values = smooth_periodic_jax(parameters, smoothing_passes)
    phase = bounded_phase_jax(values, bounds)
    if quantization_levels is not None:
        phase = quantize_phase_jax(
            phase,
            quantization_levels,
            bounds,
            straight_through=straight_through,
        )
    return phase


def _wrapped_difference_jax(values):
    return jnp.arctan2(jnp.sin(values), jnp.cos(values))


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
