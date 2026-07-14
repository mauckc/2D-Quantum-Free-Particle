"""Optional array-native JAX propagation for differentiation and JIT."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import jax
import jax.numpy as jnp


JaxArray = jax.Array

__all__ = [
    "JaxPropagationBackend",
    "apply_phase_jax",
    "fresnel_transfer_jax",
    "multiplane_scan_jax",
    "propagate_fresnel_jax",
]


@dataclass(frozen=True)
class JaxPropagationBackend:
    """JAX implementation of the generic propagation backend protocol."""

    dtype: str = "complex128"
    name: str = "jax"

    def asarray(self, values: object) -> JaxArray:
        return jnp.asarray(values, dtype=getattr(jnp, self.dtype))

    def shape(self, values: JaxArray) -> tuple[int, ...]:
        return values.shape

    def exp(self, values: JaxArray) -> JaxArray:
        return jnp.exp(values)

    def multiply(
        self,
        left: JaxArray,
        right: JaxArray | float | complex,
    ) -> JaxArray:
        return jnp.multiply(left, right)

    def fft2(self, values: JaxArray) -> JaxArray:
        return jnp.fft.fft2(values)

    def ifft2(self, values: JaxArray) -> JaxArray:
        return jnp.fft.ifft2(values)


def fresnel_transfer_jax(
    fx: Any,
    fy: Any,
    wavelength: float | JaxArray,
    distance: float | JaxArray,
) -> JaxArray:
    """Create a backend-native paraxial transfer function.

    ``fx`` and ``fy`` are broadcast-compatible spatial-frequency arrays in
    cycles per unit length. No traced value is converted through NumPy.
    """

    fx_values = jnp.asarray(fx)
    fy_values = jnp.asarray(fy)
    wavelength_value = jnp.asarray(wavelength, dtype=fx_values.real.dtype)
    distance_value = jnp.asarray(distance, dtype=fx_values.real.dtype)
    wave_number = 2.0 * jnp.pi / wavelength_value
    phase = wave_number * distance_value - jnp.pi * wavelength_value * distance_value * (
        fx_values**2 + fy_values**2
    )
    return jnp.exp(1j * phase)


def propagate_fresnel_jax(field: Any, transfer: Any) -> JaxArray:
    """Apply a prepared Fresnel transfer without host-array conversion."""

    values = jnp.asarray(field)
    transfer_values = jnp.asarray(transfer)
    return jnp.fft.ifft2(jnp.fft.fft2(values) * transfer_values)


def apply_phase_jax(field: Any, phase: Any) -> JaxArray:
    """Apply a differentiable unit-magnitude phase-only element."""

    return jnp.asarray(field) * jnp.exp(1j * jnp.asarray(phase))


def multiplane_scan_jax(
    field: Any,
    phase_masks: Any,
    transfers: Any,
) -> tuple[JaxArray, JaxArray]:
    """Apply phase/propagation pairs using a fixed-shape JAX scan.

    Returns the final field and the field after every propagation plane. All
    planes must share a transverse shape, which keeps compilation static.
    """

    initial = jnp.asarray(field)
    phases = jnp.asarray(phase_masks)
    transfer_values = jnp.asarray(transfers)

    def step(state: JaxArray, data: tuple[JaxArray, JaxArray]):
        phase, transfer = data
        next_state = propagate_fresnel_jax(apply_phase_jax(state, phase), transfer)
        return next_state, next_state

    return jax.lax.scan(step, initial, (phases, transfer_values))
