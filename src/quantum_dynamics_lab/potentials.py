from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .config import PotentialConfig
from .grid import Grid


FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def build_potential(config: PotentialConfig, grid: Grid) -> FloatArray:
    kind = config.kind.lower()
    if kind == "free":
        return np.zeros((grid.n, grid.n), dtype=np.float64)
    if kind == "barrier":
        return np.where(_barrier_mask(config, grid), config.height, 0.0).astype(np.float64)
    if kind == "double_slit":
        wall = _barrier_mask(config, grid)
        slit = upper_slit_mask(config, grid) | lower_slit_mask(config, grid)
        return np.where(wall & ~slit, config.height, 0.0).astype(np.float64)
    raise ValueError(f"unsupported potential kind: {config.kind}")


def _barrier_mask(config: PotentialConfig, grid: Grid) -> BoolArray:
    return np.abs(grid.X - config.barrier_x) <= config.barrier_width / 2.0


def transmission_mask(config: PotentialConfig, grid: Grid) -> BoolArray:
    return grid.X > config.barrier_x + config.barrier_width / 2.0


def upper_slit_mask(config: PotentialConfig, grid: Grid) -> BoolArray:
    half_width = config.slit_width / 2.0
    center = config.slit_separation / 2.0
    return _barrier_mask(config, grid) & (np.abs(grid.Y - center) <= half_width)


def lower_slit_mask(config: PotentialConfig, grid: Grid) -> BoolArray:
    half_width = config.slit_width / 2.0
    center = -config.slit_separation / 2.0
    return _barrier_mask(config, grid) & (np.abs(grid.Y - center) <= half_width)


def which_path_masks(config: PotentialConfig, grid: Grid) -> tuple[BoolArray, BoolArray]:
    return grid.Y >= 0.0, grid.Y < 0.0
