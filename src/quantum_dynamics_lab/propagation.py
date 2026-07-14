"""Pure, backend-neutral split-step propagation primitives.

The caller owns physical units, grid construction, generator preparation,
boundaries, normalization, diagnostics, and persistence. This module owns only
the symmetric position/spectral composition and delegates every array operation
to an explicit backend.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Generic, NamedTuple, Protocol, TypeVar

import numpy as np
from numpy.typing import NDArray


ArrayT = TypeVar("ArrayT")
ComplexArray = NDArray[np.complex128]
Scalar = float | complex

__all__ = [
    "NumPyPropagationBackend",
    "PropagationBackend",
    "SplitStepOperators",
    "prepare_split_step_operators",
    "propagate_split_step",
]


class PropagationBackend(Protocol[ArrayT]):
    """Array operations required by the split-step propagation kernel.

    A backend owns array conversion, complex dtype/device policy, elementwise
    operations, and FFT execution. Implementations must not convert through a
    different array system inside these methods.
    """

    name: str

    def asarray(self, values: object) -> ArrayT:
        ...

    def shape(self, values: ArrayT) -> tuple[int, ...]:
        ...

    def exp(self, values: ArrayT) -> ArrayT:
        ...

    def multiply(self, left: ArrayT, right: ArrayT | Scalar) -> ArrayT:
        ...

    def fft2(self, values: ArrayT) -> ArrayT:
        ...

    def ifft2(self, values: ArrayT) -> ArrayT:
        ...


@dataclass(frozen=True)
class NumPyPropagationBackend:
    """Complex128 NumPy correctness backend for propagation."""

    name: str = "numpy"

    def asarray(self, values: object) -> ComplexArray:
        return np.asarray(values, dtype=np.complex128)

    def shape(self, values: ComplexArray) -> tuple[int, ...]:
        return values.shape

    def exp(self, values: ComplexArray) -> ComplexArray:
        return np.asarray(np.exp(values), dtype=np.complex128)

    def multiply(
        self, left: ComplexArray, right: ComplexArray | Scalar
    ) -> ComplexArray:
        return np.asarray(np.multiply(left, right), dtype=np.complex128)

    def fft2(self, values: ComplexArray) -> ComplexArray:
        return np.asarray(np.fft.fft2(values), dtype=np.complex128)

    def ifft2(self, values: ComplexArray) -> ComplexArray:
        return np.asarray(np.fft.ifft2(values), dtype=np.complex128)


class SplitStepOperators(NamedTuple, Generic[ArrayT]):
    """Backend-native phase arrays for one symmetric split step."""

    position_half_step: ArrayT
    spectral_step: ArrayT


def prepare_split_step_operators(
    position_generator: object,
    spectral_generator: object,
    step_size: float,
    *,
    backend: PropagationBackend[ArrayT],
) -> SplitStepOperators[ArrayT]:
    """Exponentiate position and spectral generators for one step.

    Generators are rates in the evolution coordinate: the position generator is
    applied for half a step on either side of the full spectral-generator step.
    Negative finite step sizes are supported for reversibility checks.
    """

    if not math.isfinite(step_size):
        raise ValueError("step_size must be finite")
    position = backend.asarray(position_generator)
    spectral = backend.asarray(spectral_generator)
    _validate_generator_shapes(position, spectral, backend)
    return SplitStepOperators(
        position_half_step=backend.exp(
            backend.multiply(position, 0.5 * step_size)
        ),
        spectral_step=backend.exp(backend.multiply(spectral, step_size)),
    )


def propagate_split_step(
    field: object,
    operators: SplitStepOperators[ArrayT],
    *,
    backend: PropagationBackend[ArrayT],
) -> ArrayT:
    """Advance one two-dimensional field with a symmetric split step.

    The function is pure: it returns a new backend array and does not mutate the
    field or operators. The caller is responsible for applying boundaries and
    computing diagnostics after this periodic FFT propagation step.
    """

    state = backend.asarray(field)
    position = backend.asarray(operators.position_half_step)
    spectral = backend.asarray(operators.spectral_step)
    _validate_propagation_shapes(state, position, spectral, backend)

    state = backend.multiply(position, state)
    spectrum = backend.fft2(state)
    state = backend.ifft2(backend.multiply(spectrum, spectral))
    return backend.multiply(position, state)


def _validate_generator_shapes(
    position: ArrayT,
    spectral: ArrayT,
    backend: PropagationBackend[ArrayT],
) -> None:
    position_shape = backend.shape(position)
    spectral_shape = backend.shape(spectral)
    if len(position_shape) != 2 or len(spectral_shape) != 2:
        raise ValueError("split-step generators must be two-dimensional")
    if position_shape != spectral_shape:
        raise ValueError("position and spectral generators must have matching shapes")


def _validate_propagation_shapes(
    field: ArrayT,
    position: ArrayT,
    spectral: ArrayT,
    backend: PropagationBackend[ArrayT],
) -> None:
    shapes = (
        backend.shape(field),
        backend.shape(position),
        backend.shape(spectral),
    )
    if any(len(shape) != 2 for shape in shapes):
        raise ValueError("field and split-step operators must be two-dimensional")
    if len(set(shapes)) != 1:
        raise ValueError("field and split-step operators must have matching shapes")
