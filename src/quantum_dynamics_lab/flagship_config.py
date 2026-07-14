"""Configuration for the robust Gaussian-to-HG10 flagship experiment."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
import math
import tomllib

from .wave_config import OpticalGridConfig


@dataclass(frozen=True)
class PerturbationScenario:
    name: str
    kind: str
    wavelength_scale: float = 1.0
    alignment_pixels: tuple[int, int] = (0, 0)
    phase_depth_scale: float = 1.0
    spacing_scale: float = 1.0

    @property
    def signature(self) -> tuple[Any, ...]:
        return (
            self.kind,
            self.wavelength_scale,
            self.alignment_pixels,
            self.phase_depth_scale,
            self.spacing_scale,
        )


@dataclass(frozen=True)
class FlagshipSystemConfig:
    waist_m: float = 0.55e-3
    plane_distances_m: tuple[float, float, float] = (0.04, 0.05, 0.03)
    aperture_radius_m: float = 1.5e-3


@dataclass(frozen=True)
class FlagshipConstraintConfig:
    phase_bounds_rad: tuple[float, float] = (-math.pi, math.pi)
    smoothing_passes: int = 1
    total_variation_weight: float = 1.0e-4


@dataclass(frozen=True)
class FlagshipObjectiveConfig:
    intensity_weight: float = 0.05
    worst_case_weight: float = 0.5


@dataclass(frozen=True)
class FlagshipOptimizerConfig:
    iterations: int = 60
    learning_rate: float = 0.08
    beta1: float = 0.9
    beta2: float = 0.999
    epsilon: float = 1.0e-8
    initial_scale: float = 0.2


@dataclass(frozen=True)
class FlagshipConfig:
    name: str = "robust_gaussian_to_hg10"
    seed: int = 20260714
    grid: OpticalGridConfig = field(default_factory=OpticalGridConfig)
    system: FlagshipSystemConfig = field(default_factory=FlagshipSystemConfig)
    constraints: FlagshipConstraintConfig = field(
        default_factory=FlagshipConstraintConfig
    )
    objective: FlagshipObjectiveConfig = field(default_factory=FlagshipObjectiveConfig)
    optimizer: FlagshipOptimizerConfig = field(default_factory=FlagshipOptimizerConfig)
    optimization_scenarios: tuple[PerturbationScenario, ...] = ()
    held_out_scenarios: tuple[PerturbationScenario, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_flagship_config(path: str | Path) -> FlagshipConfig:
    path = Path(path)
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    grid_data = _table(data, "grid")
    system_data = _table(data, "system")
    constraint_data = _table(data, "constraints")
    objective_data = _table(data, "objective")
    optimizer_data = _table(data, "optimizer")
    robustness_data = _table(data, "robustness")
    distances = tuple(
        float(value)
        for value in system_data.get("plane_distances_m", [0.04, 0.05, 0.03])
    )
    if len(distances) != 3:
        raise ValueError("flagship requires exactly three plane distances")
    return FlagshipConfig(
        name=str(data.get("name", path.stem)),
        seed=int(data.get("seed", 20260714)),
        grid=OpticalGridConfig(
            shape=_int_pair(grid_data.get("shape", 32), "grid.shape"),
            extent_m=_float_pair(grid_data.get("extent_m", 4.0e-3), "grid.extent_m"),
            wavelength_m=float(grid_data.get("wavelength_m", 1064.0e-9)),
        ),
        system=FlagshipSystemConfig(
            waist_m=float(system_data.get("waist_m", 0.55e-3)),
            plane_distances_m=distances,
            aperture_radius_m=float(
                system_data.get("aperture_radius_m", 1.5e-3)
            ),
        ),
        constraints=FlagshipConstraintConfig(
            phase_bounds_rad=_float_pair(
                constraint_data.get("phase_bounds_rad", (-math.pi, math.pi)),
                "constraints.phase_bounds_rad",
            ),
            smoothing_passes=int(constraint_data.get("smoothing_passes", 1)),
            total_variation_weight=float(
                constraint_data.get("total_variation_weight", 1.0e-4)
            ),
        ),
        objective=FlagshipObjectiveConfig(
            intensity_weight=float(objective_data.get("intensity_weight", 0.05)),
            worst_case_weight=float(objective_data.get("worst_case_weight", 0.5)),
        ),
        optimizer=FlagshipOptimizerConfig(
            iterations=int(optimizer_data.get("iterations", 60)),
            learning_rate=float(optimizer_data.get("learning_rate", 0.08)),
            beta1=float(optimizer_data.get("beta1", 0.9)),
            beta2=float(optimizer_data.get("beta2", 0.999)),
            epsilon=float(optimizer_data.get("epsilon", 1.0e-8)),
            initial_scale=float(optimizer_data.get("initial_scale", 0.2)),
        ),
        optimization_scenarios=_scenarios(
            robustness_data.get("optimization_scenarios", [])
        ),
        held_out_scenarios=_scenarios(
            robustness_data.get("held_out_scenarios", [])
        ),
    )


def _scenarios(items: Any) -> tuple[PerturbationScenario, ...]:
    if not isinstance(items, list):
        raise ValueError("scenario sets must be arrays of inline tables")
    scenarios = []
    for item in items:
        if not isinstance(item, dict) or "name" not in item or "kind" not in item:
            raise ValueError("each scenario requires name and kind")
        scenarios.append(
            PerturbationScenario(
                name=str(item["name"]),
                kind=str(item["kind"]),
                wavelength_scale=float(item.get("wavelength_scale", 1.0)),
                alignment_pixels=_int_pair(
                    item.get("alignment_pixels", (0, 0)),
                    "alignment_pixels",
                ),
                phase_depth_scale=float(item.get("phase_depth_scale", 1.0)),
                spacing_scale=float(item.get("spacing_scale", 1.0)),
            )
        )
    return tuple(scenarios)


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
