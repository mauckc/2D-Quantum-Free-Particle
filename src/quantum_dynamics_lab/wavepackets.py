from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .config import WavePacketConfig
from .grid import Grid, normalize


def gaussian_packet(config: WavePacketConfig, grid: Grid) -> NDArray[np.complex128]:
    cx, cy = config.center
    px, py = config.momentum
    sigma = config.sigma
    if sigma <= 0:
        raise ValueError("wave_packet.sigma must be positive")
    envelope = np.exp(
        -(((grid.X - cx) ** 2) + ((grid.Y - cy) ** 2)) / (4.0 * sigma * sigma)
    )
    phase = np.exp(1j * (px * grid.X + py * grid.Y))
    psi = config.amplitude * envelope * phase
    return normalize(psi.astype(np.complex128), grid)
