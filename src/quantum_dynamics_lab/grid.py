from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .config import GridConfig


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Grid:
    n: int
    length: float
    dx: float
    dy: float
    x: FloatArray
    y: FloatArray
    X: FloatArray
    Y: FloatArray
    kx: FloatArray
    ky: FloatArray
    k2: FloatArray


def make_grid(config: GridConfig) -> Grid:
    if config.n < 4:
        raise ValueError("grid.n must be at least 4")
    if config.length <= 0:
        raise ValueError("grid.length must be positive")
    dx = config.length / config.n
    x = (np.arange(config.n, dtype=np.float64) - config.n / 2) * dx
    y = x.copy()
    X, Y = np.meshgrid(x, y, indexing="xy")
    k = 2.0 * np.pi * np.fft.fftfreq(config.n, d=dx)
    kx, ky = np.meshgrid(k, k, indexing="xy")
    return Grid(
        n=config.n,
        length=config.length,
        dx=dx,
        dy=dx,
        x=x,
        y=y,
        X=X,
        Y=Y,
        kx=kx,
        ky=ky,
        k2=kx * kx + ky * ky,
    )


def integrate_probability(probability: FloatArray, grid: Grid) -> float:
    return float(np.sum(probability) * grid.dx * grid.dy)


def norm(psi: NDArray[np.complex128], grid: Grid) -> float:
    return integrate_probability(np.abs(psi) ** 2, grid)


def normalize(psi: NDArray[np.complex128], grid: Grid) -> NDArray[np.complex128]:
    current = norm(psi, grid)
    if current <= 0:
        raise ValueError("cannot normalize a zero-norm wavefunction")
    return psi / np.sqrt(current)
