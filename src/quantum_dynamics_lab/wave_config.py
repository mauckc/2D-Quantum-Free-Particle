"""Typed TOML configuration for scalar-optics forward runs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
import tomllib


@dataclass(frozen=True)
class OpticalGridConfig:
    shape: tuple[int, int] = (128, 128)
    extent_m: tuple[float, float] = (6.0e-3, 6.0e-3)
    wavelength_m: float = 1064.0e-9


@dataclass(frozen=True)
class OpticalSourceConfig:
    kind: str = "gaussian"
    waist_m: float = 0.5e-3
    center_m: tuple[float, float] = (0.0, 0.0)
    order: tuple[int, int] = (0, 0)
    path: Path | None = None
    array_key: str = "field"


@dataclass(frozen=True)
class OpticalStepConfig:
    kind: str
    model: str = "fresnel"
    distance_m: float | None = None
    focal_length_m: float | None = None
    radius_m: float | None = None
    size_m: tuple[float, float] | None = None
    center_m: tuple[float, float] = (0.0, 0.0)
    phase_profile: str = "sinusoidal"
    phase_amplitude_rad: float = 0.0
    phase_period_m: float | None = None
    axis: str = "x"


@dataclass(frozen=True)
class WavePropagationConfig:
    name: str = "wave_propagation"
    grid: OpticalGridConfig = field(default_factory=OpticalGridConfig)
    source: OpticalSourceConfig = field(default_factory=OpticalSourceConfig)
    steps: tuple[OpticalStepConfig, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        if self.source.path is not None:
            data["source"]["path"] = str(self.source.path)
        return data


def load_wave_config(path: str | Path) -> WavePropagationConfig:
    """Load an optics run without sharing the legacy quantum config schema."""

    path = Path(path).resolve()
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    grid_data = _table(data, "grid")
    source_data = _table(data, "source")
    source_path = source_data.get("path")
    resolved_source_path = None
    if source_path is not None:
        resolved_source_path = Path(str(source_path))
        if not resolved_source_path.is_absolute():
            resolved_source_path = path.parent / resolved_source_path
        resolved_source_path = resolved_source_path.resolve()

    raw_steps = data.get("steps", [])
    if not isinstance(raw_steps, list):
        raise ValueError("[[steps]] must be an array of TOML tables")
    steps = []
    for index, item in enumerate(raw_steps):
        if not isinstance(item, dict) or "kind" not in item:
            raise ValueError(f"steps[{index}] must be a table with a kind")
        raw_size = item.get("size_m")
        steps.append(
            OpticalStepConfig(
                kind=str(item["kind"]),
                model=str(item.get("model", "fresnel")),
                distance_m=_optional_float(item.get("distance_m")),
                focal_length_m=_optional_float(item.get("focal_length_m")),
                radius_m=_optional_float(item.get("radius_m")),
                size_m=None if raw_size is None else _float_pair(raw_size, "size_m"),
                center_m=_float_pair(
                    item.get("center_m", (0.0, 0.0)),
                    "center_m",
                ),
                phase_profile=str(item.get("phase_profile", "sinusoidal")),
                phase_amplitude_rad=float(item.get("phase_amplitude_rad", 0.0)),
                phase_period_m=_optional_float(item.get("phase_period_m")),
                axis=str(item.get("axis", "x")),
            )
        )

    return WavePropagationConfig(
        name=str(data.get("name", path.stem)),
        grid=OpticalGridConfig(
            shape=_int_pair(grid_data.get("shape", 128), "grid.shape"),
            extent_m=_float_pair(
                grid_data.get("extent_m", 6.0e-3),
                "grid.extent_m",
            ),
            wavelength_m=float(grid_data.get("wavelength_m", 1064.0e-9)),
        ),
        source=OpticalSourceConfig(
            kind=str(source_data.get("kind", "gaussian")),
            waist_m=float(source_data.get("waist_m", 0.5e-3)),
            center_m=_float_pair(
                source_data.get("center_m", (0.0, 0.0)),
                "source.center_m",
            ),
            order=_int_pair(source_data.get("order", (0, 0)), "source.order"),
            path=resolved_source_path,
            array_key=str(source_data.get("array_key", "field")),
        ),
        steps=tuple(steps),
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


def _optional_float(value: Any) -> float | None:
    return None if value is None else float(value)
