"""Typed configuration for constrained scalar-wave inverse design."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
import math
import tomllib

from .wave_config import OpticalGridConfig


@dataclass(frozen=True)
class DesignModeConfig:
    kind: str = "gaussian"
    waist_m: float = 0.5e-3
    order: tuple[int, int] = (0, 0)
    center_m: tuple[float, float] = (0.0, 0.0)


@dataclass(frozen=True)
class PhaseDesignConfig:
    plane_distances_m: tuple[float, ...] = (0.05,)
    initial_scale: float = 0.1
    phase_bounds_rad: tuple[float, float] = (-math.pi, math.pi)
    smoothing_passes: int = 1
    total_variation_weight: float = 0.0
    quantization_levels: int | None = None


@dataclass(frozen=True)
class ObjectiveConfig:
    mode_overlap_weight: float = 1.0
    intensity_weight: float = 0.0


@dataclass(frozen=True)
class AdamConfig:
    iterations: int = 40
    learning_rate: float = 0.05
    beta1: float = 0.9
    beta2: float = 0.999
    epsilon: float = 1.0e-8
    checkpoint_interval: int = 10
    early_stopping_patience: int = 20
    min_delta: float = 1.0e-8


@dataclass(frozen=True)
class InverseDesignConfig:
    name: str = "inverse_design"
    seed: int = 0
    grid: OpticalGridConfig = field(default_factory=OpticalGridConfig)
    input_mode: DesignModeConfig = field(default_factory=DesignModeConfig)
    target_mode: DesignModeConfig = field(
        default_factory=lambda: DesignModeConfig(kind="hermite_gaussian", order=(1, 0))
    )
    design: PhaseDesignConfig = field(default_factory=PhaseDesignConfig)
    objective: ObjectiveConfig = field(default_factory=ObjectiveConfig)
    optimizer: AdamConfig = field(default_factory=AdamConfig)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_inverse_design_config(path: str | Path) -> InverseDesignConfig:
    path = Path(path)
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    grid_data = _table(data, "grid")
    input_data = _table(data, "input_mode")
    target_data = _table(data, "target_mode")
    design_data = _table(data, "design")
    objective_data = _table(data, "objective")
    optimizer_data = _table(data, "optimizer")
    quantization_levels = design_data.get("quantization_levels")
    return InverseDesignConfig(
        name=str(data.get("name", path.stem)),
        seed=int(data.get("seed", 0)),
        grid=OpticalGridConfig(
            shape=_int_pair(grid_data.get("shape", 128), "grid.shape"),
            extent_m=_float_pair(
                grid_data.get("extent_m", 6.0e-3),
                "grid.extent_m",
            ),
            wavelength_m=float(grid_data.get("wavelength_m", 1064.0e-9)),
        ),
        input_mode=_mode(input_data, "gaussian"),
        target_mode=_mode(target_data, "hermite_gaussian", (1, 0)),
        design=PhaseDesignConfig(
            plane_distances_m=tuple(
                float(value)
                for value in design_data.get("plane_distances_m", [0.05])
            ),
            initial_scale=float(design_data.get("initial_scale", 0.1)),
            phase_bounds_rad=_float_pair(
                design_data.get("phase_bounds_rad", (-math.pi, math.pi)),
                "design.phase_bounds_rad",
            ),
            smoothing_passes=int(design_data.get("smoothing_passes", 1)),
            total_variation_weight=float(
                design_data.get("total_variation_weight", 0.0)
            ),
            quantization_levels=(
                None if quantization_levels is None else int(quantization_levels)
            ),
        ),
        objective=ObjectiveConfig(
            mode_overlap_weight=float(
                objective_data.get("mode_overlap_weight", 1.0)
            ),
            intensity_weight=float(objective_data.get("intensity_weight", 0.0)),
        ),
        optimizer=AdamConfig(
            iterations=int(optimizer_data.get("iterations", 40)),
            learning_rate=float(optimizer_data.get("learning_rate", 0.05)),
            beta1=float(optimizer_data.get("beta1", 0.9)),
            beta2=float(optimizer_data.get("beta2", 0.999)),
            epsilon=float(optimizer_data.get("epsilon", 1.0e-8)),
            checkpoint_interval=int(
                optimizer_data.get("checkpoint_interval", 10)
            ),
            early_stopping_patience=int(
                optimizer_data.get("early_stopping_patience", 20)
            ),
            min_delta=float(optimizer_data.get("min_delta", 1.0e-8)),
        ),
    )


def _mode(
    data: dict[str, Any],
    default_kind: str,
    default_order: tuple[int, int] = (0, 0),
) -> DesignModeConfig:
    return DesignModeConfig(
        kind=str(data.get("kind", default_kind)),
        waist_m=float(data.get("waist_m", 0.5e-3)),
        order=_int_pair(data.get("order", default_order), "mode.order"),
        center_m=_float_pair(data.get("center_m", (0.0, 0.0)), "mode.center_m"),
    )


def _table(data: dict[str, Any], name: str) -> dict[str, Any]:
    value = data.get(name, {})
    if not isinstance(value, dict):
        raise ValueError(f"[{name}] must be a TOML table")
    return value


def _int_pair(value: Any, name: str) -> tuple[int, int]:
    values = (value, value) if isinstance(value, int) else value
    if not isinstance(values, (list, tuple)) or len(values) != 2:
        raise ValueError(f"{name} must be an integer or two-item array")
    return int(values[0]), int(values[1])


def _float_pair(value: Any, name: str) -> tuple[float, float]:
    values = (value, value) if isinstance(value, (int, float)) else value
    if not isinstance(values, (list, tuple)) or len(values) != 2:
        raise ValueError(f"{name} must be a number or two-item array")
    return float(values[0]), float(values[1])
