from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any
import tomllib


@dataclass(frozen=True)
class GridConfig:
    n: int = 128
    length: float = 12.0


@dataclass(frozen=True)
class WavePacketConfig:
    center: tuple[float, float] = (-3.0, 0.0)
    sigma: float = 0.45
    momentum: tuple[float, float] = (5.0, 0.0)
    amplitude: float = 1.0


@dataclass(frozen=True)
class PotentialConfig:
    kind: str = "free"
    height: float = 0.0
    barrier_x: float = 0.0
    barrier_width: float = 0.25
    slit_width: float = 0.45
    slit_separation: float = 1.2


@dataclass(frozen=True)
class SolverConfig:
    dt: float = 0.002
    tf: float = 1.0
    frame_interval: float = 0.02
    hbar: float = 1.0
    mass: float = 1.0
    backend: str = "auto"


@dataclass(frozen=True)
class BoundaryConfig:
    kind: str = "periodic"
    width: float = 0.0
    strength: float = 0.0
    power: float = 2.0


@dataclass(frozen=True)
class MeasurementConfig:
    kind: str = "none"
    interval: float | None = None
    time: float | None = None
    detector_x: float | None = None
    detector_buffer: float = 0.0


@dataclass(frozen=True)
class OutputConfig:
    fps: int = 20
    render_video: bool = True


@dataclass(frozen=True)
class ExperimentConfig:
    name: str = "free_packet"
    experiment: str = "free_packet"
    grid: GridConfig = field(default_factory=GridConfig)
    wave_packet: WavePacketConfig = field(default_factory=WavePacketConfig)
    potential: PotentialConfig = field(default_factory=PotentialConfig)
    solver: SolverConfig = field(default_factory=SolverConfig)
    boundary: BoundaryConfig = field(default_factory=BoundaryConfig)
    measurement: MeasurementConfig = field(default_factory=MeasurementConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SweepVariant:
    name: str
    measurement_kind: str = "zeno"
    measurement_interval: float | None = None
    potential_height: float | None = None


@dataclass(frozen=True)
class SweepConfig:
    name: str
    base_config: Path
    variants: tuple[SweepVariant, ...]


def _tuple2(value: Any, field_name: str) -> tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{field_name} must be a two-item array")
    return (float(value[0]), float(value[1]))


def _section(data: dict[str, Any], name: str) -> dict[str, Any]:
    value = data.get(name, {})
    if not isinstance(value, dict):
        raise ValueError(f"[{name}] must be a TOML table")
    return value


def load_experiment_config(path: str | Path) -> ExperimentConfig:
    path = Path(path)
    data = tomllib.loads(path.read_text(encoding="utf-8"))

    grid_data = _section(data, "grid")
    wave_data = _section(data, "wave_packet")
    potential_data = _section(data, "potential")
    solver_data = _section(data, "solver")
    boundary_data = _section(data, "boundary")
    measurement_data = _section(data, "measurement")
    output_data = _section(data, "output")

    grid = GridConfig(
        n=int(grid_data.get("n", GridConfig.n)),
        length=float(grid_data.get("length", GridConfig.length)),
    )
    wave = WavePacketConfig(
        center=_tuple2(wave_data.get("center", WavePacketConfig.center), "wave_packet.center"),
        sigma=float(wave_data.get("sigma", WavePacketConfig.sigma)),
        momentum=_tuple2(
            wave_data.get("momentum", WavePacketConfig.momentum),
            "wave_packet.momentum",
        ),
        amplitude=float(wave_data.get("amplitude", WavePacketConfig.amplitude)),
    )
    potential = PotentialConfig(
        kind=str(potential_data.get("kind", PotentialConfig.kind)),
        height=float(potential_data.get("height", PotentialConfig.height)),
        barrier_x=float(potential_data.get("barrier_x", PotentialConfig.barrier_x)),
        barrier_width=float(potential_data.get("barrier_width", PotentialConfig.barrier_width)),
        slit_width=float(potential_data.get("slit_width", PotentialConfig.slit_width)),
        slit_separation=float(
            potential_data.get("slit_separation", PotentialConfig.slit_separation)
        ),
    )
    solver = SolverConfig(
        dt=float(solver_data.get("dt", SolverConfig.dt)),
        tf=float(solver_data.get("tf", SolverConfig.tf)),
        frame_interval=float(
            solver_data.get("frame_interval", SolverConfig.frame_interval)
        ),
        hbar=float(solver_data.get("hbar", SolverConfig.hbar)),
        mass=float(solver_data.get("mass", SolverConfig.mass)),
        backend=str(solver_data.get("backend", SolverConfig.backend)),
    )
    boundary = BoundaryConfig(
        kind=str(boundary_data.get("kind", BoundaryConfig.kind)),
        width=float(boundary_data.get("width", BoundaryConfig.width)),
        strength=float(boundary_data.get("strength", BoundaryConfig.strength)),
        power=float(boundary_data.get("power", BoundaryConfig.power)),
    )
    measurement = MeasurementConfig(
        kind=str(measurement_data.get("kind", MeasurementConfig.kind)),
        interval=(
            None
            if measurement_data.get("interval") is None
            else float(measurement_data["interval"])
        ),
        time=(
            None if measurement_data.get("time") is None else float(measurement_data["time"])
        ),
        detector_x=(
            None
            if measurement_data.get("detector_x") is None
            else float(measurement_data["detector_x"])
        ),
        detector_buffer=float(
            measurement_data.get("detector_buffer", MeasurementConfig.detector_buffer)
        ),
    )
    output = OutputConfig(
        fps=int(output_data.get("fps", OutputConfig.fps)),
        render_video=bool(output_data.get("render_video", OutputConfig.render_video)),
    )
    return ExperimentConfig(
        name=str(data.get("name", path.stem)),
        experiment=str(data.get("experiment", "free_packet")),
        grid=grid,
        wave_packet=wave,
        potential=potential,
        solver=solver,
        boundary=boundary,
        measurement=measurement,
        output=output,
    )


def load_sweep_config(path: str | Path) -> SweepConfig:
    path = Path(path)
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    variants = []
    for item in data.get("variants", []):
        interval = item.get("measurement_interval")
        height = item.get("potential_height")
        variants.append(
            SweepVariant(
                name=str(item["name"]),
                measurement_kind=str(item.get("measurement_kind", "zeno")),
                measurement_interval=None if interval is None else float(interval),
                potential_height=None if height is None else float(height),
            )
        )
    scan = data.get("scan", {})
    if scan:
        intervals = [float(value) for value in scan.get("measurement_intervals", [])]
        heights = [float(value) for value in scan.get("potential_heights", [])]
        if not intervals:
            intervals = [0.0]
        if not heights:
            heights = [0.0]
        for height in heights:
            for interval in intervals:
                variants.append(
                    SweepVariant(
                        name=_scan_variant_name(height, interval),
                        measurement_kind="none" if interval <= 0 else "zeno",
                        measurement_interval=None if interval <= 0 else interval,
                        potential_height=height,
                    )
                )
    base_config = Path(str(data["base_config"]))
    if not base_config.is_absolute():
        base_config = path.parent / base_config
    return SweepConfig(
        name=str(data.get("name", path.stem)),
        base_config=base_config,
        variants=tuple(variants),
    )


def apply_sweep_variant(config: ExperimentConfig, variant: SweepVariant) -> ExperimentConfig:
    kind = variant.measurement_kind
    interval = variant.measurement_interval
    if interval is None or interval <= 0:
        kind = "none"
        interval = None
    measurement = replace(config.measurement, kind=kind, interval=interval)
    potential = config.potential
    if variant.potential_height is not None:
        potential = replace(potential, height=variant.potential_height)
    return replace(config, name=variant.name, potential=potential, measurement=measurement)


def _scan_variant_name(height: float, interval: float) -> str:
    height_label = _number_label(height)
    if interval <= 0:
        return f"h{height_label}_unmeasured"
    return f"h{height_label}_dt{_number_label(interval)}"


def _number_label(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")
