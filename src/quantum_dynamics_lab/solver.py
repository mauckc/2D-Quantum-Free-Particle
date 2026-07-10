from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .config import SolverConfig
from .grid import Grid


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class SplitStepSolver:
    grid: Grid
    potential: FloatArray
    config: SolverConfig

    def __post_init__(self) -> None:
        if self.config.dt <= 0:
            raise ValueError("solver.dt must be positive")
        if self.config.tf <= 0:
            raise ValueError("solver.tf must be positive")
        if self.config.mass <= 0:
            raise ValueError("solver.mass must be positive")
        if self.config.hbar <= 0:
            raise ValueError("solver.hbar must be positive")

    @property
    def step_count(self) -> int:
        return max(1, int(round(self.config.tf / self.config.dt)))

    @property
    def frame_stride(self) -> int:
        return max(1, int(round(self.config.frame_interval / self.config.dt)))

    def step(self, psi: ComplexArray) -> ComplexArray:
        dt = self.config.dt
        hbar = self.config.hbar
        mass = self.config.mass
        position_half = np.exp(-0.5j * self.potential * dt / hbar)
        kinetic = np.exp(-1j * hbar * self.grid.k2 * dt / (2.0 * mass))
        psi = position_half * psi
        psi = np.fft.ifft2(np.fft.fft2(psi) * kinetic)
        return (position_half * psi).astype(np.complex128)
