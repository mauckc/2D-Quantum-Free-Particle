"""Optional JAX-differentiable scalar optical objectives."""

from __future__ import annotations

from typing import Any

import jax.numpy as jnp


def mode_overlap_jax(
    field: Any,
    target: Any,
    sample_area: float = 1.0,
):
    """JAX form of normalized phase-insensitive complex-mode overlap."""

    values = jnp.asarray(field)
    target_values = jnp.asarray(target)
    area = jnp.asarray(sample_area, dtype=values.real.dtype)
    inner = jnp.vdot(target_values, values) * area
    field_power = jnp.sum(jnp.abs(values) ** 2) * area
    target_power = jnp.sum(jnp.abs(target_values) ** 2) * area
    denominator = field_power * target_power
    return jnp.where(denominator > 0.0, jnp.abs(inner) ** 2 / denominator, jnp.nan)


def intensity_similarity_jax(
    field: Any,
    target: Any,
    sample_area: float = 1.0,
):
    """JAX form of normalized intensity cosine similarity."""

    values = jnp.asarray(field)
    target_values = jnp.asarray(target)
    area = jnp.asarray(sample_area, dtype=values.real.dtype)
    intensity = jnp.abs(values) ** 2
    target_intensity = jnp.abs(target_values) ** 2
    numerator = jnp.sum(intensity * target_intensity) * area
    denominator = jnp.sqrt(
        jnp.sum(intensity**2) * area * jnp.sum(target_intensity**2) * area
    )
    return jnp.where(denominator > 0.0, numerator / denominator, jnp.nan)


def mode_power_efficiency_jax(
    field: Any,
    target: Any,
    input_power: float,
    sample_area: float = 1.0,
):
    """JAX form of input-referenced target-mode coupling efficiency."""

    values = jnp.asarray(field)
    target_values = jnp.asarray(target)
    area = jnp.asarray(sample_area, dtype=values.real.dtype)
    input_power_value = jnp.asarray(input_power, dtype=values.real.dtype)
    target_power = jnp.sum(jnp.abs(target_values) ** 2) * area
    inner = jnp.vdot(target_values, values) * area
    denominator = target_power * input_power_value
    return jnp.where(denominator > 0.0, jnp.abs(inner) ** 2 / denominator, jnp.nan)
