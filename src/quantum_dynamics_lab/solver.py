from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import NDArray

from .backends import FFTBackend, create_backend
from .config import BoundaryConfig, SolverConfig
from .grid import Grid, norm
from .propagation import (
    PropagationBackend,
    SplitStepOperators,
    prepare_split_step_operators,
    propagate_split_step,
)


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class _QuantumPropagationBackend:
    """Adapt a legacy FFT backend to the array-native propagation protocol."""

    fft_backend: FFTBackend

    @property
    def name(self) -> str:
        return self.fft_backend.name

    def asarray(self, values: object) -> ComplexArray:
        return np.asarray(values, dtype=np.complex128)

    def shape(self, values: ComplexArray) -> tuple[int, ...]:
        return values.shape

    def exp(self, values: ComplexArray) -> ComplexArray:
        return np.asarray(np.exp(values), dtype=np.complex128)

    def multiply(
        self, left: ComplexArray, right: ComplexArray | float | complex
    ) -> ComplexArray:
        return np.asarray(np.multiply(left, right), dtype=np.complex128)

    def fft2(self, values: ComplexArray) -> ComplexArray:
        return np.asarray(self.fft_backend.fft2(values), dtype=np.complex128)

    def ifft2(self, values: ComplexArray) -> ComplexArray:
        return np.asarray(self.fft_backend.ifft2(values), dtype=np.complex128)


@dataclass(frozen=True)
class StepResult:
    psi: ComplexArray
    absorbed_probability: float


@dataclass(frozen=True)
class SplitStepSolver:
    grid: Grid
    potential: FloatArray
    config: SolverConfig
    boundary: BoundaryConfig = BoundaryConfig()
    backend: FFTBackend | None = None

    def __post_init__(self) -> None:
        if self.config.dt <= 0:
            raise ValueError("solver.dt must be positive")
        if self.config.tf <= 0:
            raise ValueError("solver.tf must be positive")
        if self.config.mass <= 0:
            raise ValueError("solver.mass must be positive")
        if self.config.hbar <= 0:
            raise ValueError("solver.hbar must be positive")
        if self.boundary.width < 0:
            raise ValueError("boundary.width must be non-negative")
        if self.boundary.strength < 0:
            raise ValueError("boundary.strength must be non-negative")
        if self.boundary.power <= 0:
            raise ValueError("boundary.power must be positive")

    @property
    def step_count(self) -> int:
        return max(1, int(round(self.config.tf / self.config.dt)))

    @property
    def frame_stride(self) -> int:
        return max(1, int(round(self.config.frame_interval / self.config.dt)))

    @cached_property
    def fft_backend(self) -> FFTBackend:
        if self.backend is not None:
            return self.backend
        return create_backend(self.config.backend)

    @property
    def backend_name(self) -> str:
        return self.fft_backend.name

    @cached_property
    def propagation_backend(self) -> PropagationBackend[ComplexArray]:
        return _QuantumPropagationBackend(self.fft_backend)

    @cached_property
    def split_step_operators(self) -> SplitStepOperators[ComplexArray]:
        return prepare_split_step_operators(
            -1j * self.potential / self.config.hbar,
            -1j * self.config.hbar * self.grid.k2 / (2.0 * self.config.mass),
            self.config.dt,
            backend=self.propagation_backend,
        )

    def step(self, psi: ComplexArray) -> ComplexArray:
        return self.step_with_diagnostics(psi).psi

    def step_with_diagnostics(self, psi: ComplexArray) -> StepResult:
        psi = propagate_split_step(
            psi,
            self.split_step_operators,
            backend=self.propagation_backend,
        )
        before_absorber = norm(psi, self.grid)
        psi = self.absorber_mask * psi
        absorbed = max(0.0, before_absorber - norm(psi, self.grid))
        return StepResult(psi=psi.astype(np.complex128), absorbed_probability=absorbed)

    @property
    def absorber_mask(self) -> FloatArray:
        if self.boundary.kind.lower() != "absorbing":
            return np.ones((self.grid.n, self.grid.n), dtype=np.float64)
        if self.boundary.width <= 0 or self.boundary.strength <= 0:
            return np.ones((self.grid.n, self.grid.n), dtype=np.float64)
        half_length = self.grid.length / 2.0
        distance_to_edge = np.minimum.reduce(
            [
                self.grid.X + half_length,
                half_length - self.grid.X,
                self.grid.Y + half_length,
                half_length - self.grid.Y,
            ]
        )
        scaled = np.clip((self.boundary.width - distance_to_edge) / self.boundary.width, 0.0, 1.0)
        damping = self.boundary.strength * np.power(scaled, self.boundary.power)
        return np.exp(-damping * self.config.dt).astype(np.float64)
