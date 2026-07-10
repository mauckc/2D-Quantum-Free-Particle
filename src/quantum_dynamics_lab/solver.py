from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .config import BoundaryConfig, SolverConfig
from .grid import Grid, norm


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


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

    def step(self, psi: ComplexArray) -> ComplexArray:
        return self.step_with_diagnostics(psi).psi

    def step_with_diagnostics(self, psi: ComplexArray) -> StepResult:
        dt = self.config.dt
        hbar = self.config.hbar
        mass = self.config.mass
        position_half = np.exp(-0.5j * self.potential * dt / hbar)
        kinetic = np.exp(-1j * hbar * self.grid.k2 * dt / (2.0 * mass))
        psi = position_half * psi
        psi = np.fft.ifft2(np.fft.fft2(psi) * kinetic)
        psi = (position_half * psi).astype(np.complex128)
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
