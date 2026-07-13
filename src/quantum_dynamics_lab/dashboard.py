from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

from .artifacts import save_run
from .config import (
    BoundaryConfig,
    ExperimentConfig,
    GridConfig,
    MeasurementConfig,
    OutputConfig,
    PotentialConfig,
    SolverConfig,
    WavePacketConfig,
)
from .experiments import ExperimentResult, run_experiment
from .render import render_run
from .validation import load_validation_config, run_validation


@dataclass(frozen=True)
class DashboardParameters:
    name: str = "barrier_zeno"
    experiment: str = "barrier_zeno"
    grid_n: int = 48
    grid_length: float = 10.0
    center_x: float = -2.5
    center_y: float = 0.0
    sigma: float = 0.5
    momentum_x: float = 4.5
    momentum_y: float = 0.0
    dt: float = 0.004
    tf: float = 0.8
    frame_interval: float = 0.04
    backend: str = "numpy"
    barrier_height: float = 40.0
    barrier_x: float = 0.0
    barrier_width: float = 0.25
    measurement_kind: str = "zeno"
    measurement_interval: float | None = 0.12
    detector_buffer: float = 0.25
    absorbing_boundary: bool = True
    boundary_width: float = 1.0
    boundary_strength: float = 4.0
    boundary_power: float = 2.0


@dataclass(frozen=True)
class DashboardRun:
    result: ExperimentResult
    run_path: Path
    outputs: dict[str, Path]


def build_experiment_config(parameters: DashboardParameters) -> ExperimentConfig:
    _validate_parameters(parameters)
    is_barrier = parameters.experiment == "barrier_zeno"
    measurement_kind = parameters.measurement_kind if is_barrier else "none"
    interval = parameters.measurement_interval if measurement_kind == "zeno" else None
    return ExperimentConfig(
        name=parameters.name,
        experiment=parameters.experiment,
        grid=GridConfig(n=parameters.grid_n, length=parameters.grid_length),
        wave_packet=WavePacketConfig(
            center=(parameters.center_x, parameters.center_y),
            sigma=parameters.sigma,
            momentum=(parameters.momentum_x, parameters.momentum_y),
        ),
        potential=PotentialConfig(
            kind="barrier" if is_barrier else "free",
            height=parameters.barrier_height if is_barrier else 0.0,
            barrier_x=parameters.barrier_x,
            barrier_width=parameters.barrier_width,
        ),
        solver=SolverConfig(
            dt=parameters.dt,
            tf=parameters.tf,
            frame_interval=parameters.frame_interval,
            backend=parameters.backend,
        ),
        boundary=BoundaryConfig(
            kind="absorbing" if parameters.absorbing_boundary else "periodic",
            width=parameters.boundary_width if parameters.absorbing_boundary else 0.0,
            strength=parameters.boundary_strength if parameters.absorbing_boundary else 0.0,
            power=parameters.boundary_power,
        ),
        measurement=MeasurementConfig(
            kind=measurement_kind,
            interval=interval,
            detector_buffer=parameters.detector_buffer if is_barrier else 0.0,
        ),
        output=OutputConfig(render_video=False),
    )


def run_dashboard_experiment(
    parameters: DashboardParameters, out_dir: str | Path
) -> DashboardRun:
    config = build_experiment_config(parameters)
    result = run_experiment(config)
    run_path = save_run(result, out_dir)
    outputs = render_run(run_path, out_dir)
    return DashboardRun(result=result, run_path=run_path, outputs=outputs)


def run_dashboard_validation(
    config_path: str | Path, out_dir: str | Path
) -> dict[str, Any]:
    return run_validation(load_validation_config(config_path), out_dir)


def timestamped_output_dir(
    root: str | Path, name: str, when: datetime | None = None
) -> Path:
    timestamp = (when or datetime.now(UTC)).strftime("%Y%m%dT%H%M%S%fZ")
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-") or "run"
    return Path(root) / f"{slug}-{timestamp}"


def experiment_summary(result: ExperimentResult) -> dict[str, float]:
    return {
        "final_norm": float(result.metrics["final_norm"]),
        "max_norm_drift": float(np.max(np.abs(result.arrays["norm"] - 1.0))),
        "final_unconditional_transmission_probability": float(
            result.metrics["final_unconditional_transmission_probability"]
        ),
        "final_detected_probability": float(result.metrics["final_detected_probability"]),
        "final_survival_weight": float(result.metrics["final_survival_weight"]),
        "final_absorbed_probability": float(result.metrics["final_absorbed_probability"]),
    }


def plot_field_panels(result: ExperimentResult, frame_index: int) -> tuple[Figure, ...]:
    frames = result.arrays["psi_frames"]
    if not 0 <= frame_index < len(frames):
        raise IndexError("frame_index is outside the saved frame range")
    psi = frames[frame_index]
    probability = np.abs(psi) ** 2
    peak = float(np.max(probability))
    phase = np.ma.masked_where(probability < peak * 1e-3, np.angle(psi))
    extent = _extent(result)

    probability_figure, probability_axis = plt.subplots(figsize=(4.4, 3.8))
    probability_image = probability_axis.imshow(
        probability, origin="lower", extent=extent, cmap="viridis", aspect="equal"
    )
    probability_axis.set_title("Probability density")
    probability_figure.colorbar(probability_image, ax=probability_axis, label="|psi|^2")

    phase_figure, phase_axis = plt.subplots(figsize=(4.4, 3.8))
    phase_image = phase_axis.imshow(
        phase,
        origin="lower",
        extent=extent,
        cmap="twilight",
        vmin=-np.pi,
        vmax=np.pi,
        aspect="equal",
    )
    phase_axis.set_title("Phase in occupied region")
    phase_figure.colorbar(phase_image, ax=phase_axis, label="radians")

    potential_figure, potential_axis = plt.subplots(figsize=(4.4, 3.8))
    potential_image = potential_axis.imshow(
        result.potential, origin="lower", extent=extent, cmap="magma", aspect="equal"
    )
    potential_axis.set_title("Potential")
    potential_figure.colorbar(potential_image, ax=potential_axis, label="V(x, y)")

    if result.config.experiment == "barrier_zeno":
        detector_x = float(result.metrics["detector_x"])
        probability_axis.axvline(
            detector_x, color="white", linestyle="--", linewidth=1.2
        )
        probability_axis.text(
            detector_x,
            extent[3],
            " detector",
            color="white",
            ha="left",
            va="top",
        )
    figures = (probability_figure, phase_figure, potential_figure)
    for figure, axis in zip(
        figures, (probability_axis, phase_axis, potential_axis), strict=True
    ):
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        figure.tight_layout()
    return figures


def plot_diagnostic_panels(result: ExperimentResult) -> tuple[Figure, ...]:
    arrays = result.arrays
    metric_times = arrays["metric_times"]
    norm_figure, norm_axis = plt.subplots(figsize=(5.4, 3.6))
    probability_figure, probability_axis = plt.subplots(figsize=(5.4, 3.6))

    norm_drift = arrays["norm"] - 1.0
    norm_axis.plot(
        metric_times, norm_drift, color="#2563eb", label="conditional norm - 1"
    )
    norm_axis.axhline(0.0, color="#6b7280", linewidth=0.8)
    norm_axis.set_title("Norm drift")
    norm_axis.set_xlabel("time")
    norm_axis.set_ylabel("probability mass")
    norm_axis.grid(True, alpha=0.25)

    if result.config.experiment == "barrier_zeno":
        series = (
            ("unconditional_transmission_probability", "no-click transmission"),
            ("detected_probability", "detector clicks"),
            ("survival_weight", "survival"),
            ("absorbed_probability", "absorbed"),
        )
    else:
        series = (
            ("survival_weight", "survival"),
            ("absorbed_probability", "absorbed"),
        )
    for key, label in series:
        probability_axis.plot(metric_times, arrays[key], label=label)
    probability_axis.set_title("Probability accounting")
    probability_axis.set_xlabel("time")
    probability_axis.set_ylabel("probability")
    probability_axis.set_ylim(-0.02, 1.02)
    probability_axis.grid(True, alpha=0.25)
    probability_axis.legend(loc="best")
    norm_figure.tight_layout()
    probability_figure.tight_layout()
    return norm_figure, probability_figure


def _validate_parameters(parameters: DashboardParameters) -> None:
    if parameters.experiment not in {"free_packet", "barrier_zeno"}:
        raise ValueError("experiment must be free_packet or barrier_zeno")
    if parameters.measurement_kind not in {"none", "zeno"}:
        raise ValueError("measurement_kind must be none or zeno")
    if parameters.grid_n < 16:
        raise ValueError("grid_n must be at least 16")
    if parameters.grid_length <= 0:
        raise ValueError("grid_length must be positive")
    half_length = parameters.grid_length / 2.0
    if not (-half_length < parameters.center_x < half_length):
        raise ValueError("packet center_x must lie inside the grid")
    if not (-half_length < parameters.center_y < half_length):
        raise ValueError("packet center_y must lie inside the grid")
    if parameters.sigma <= 0:
        raise ValueError("sigma must be positive")
    if parameters.dt <= 0 or parameters.tf <= 0:
        raise ValueError("dt and tf must be positive")
    if parameters.dt > parameters.tf:
        raise ValueError("dt must not exceed tf")
    if parameters.frame_interval < parameters.dt:
        raise ValueError("frame_interval must be at least dt")
    if parameters.backend not in {"numpy", "auto", "pyfftw", "cupy"}:
        raise ValueError("unsupported solver backend")
    if parameters.absorbing_boundary:
        if not 0 < parameters.boundary_width < half_length:
            raise ValueError("boundary_width must fit inside half the grid")
        if parameters.boundary_strength <= 0 or parameters.boundary_power <= 0:
            raise ValueError("absorbing boundary strength and power must be positive")
    if parameters.experiment != "barrier_zeno":
        return
    if parameters.barrier_width <= 0:
        raise ValueError("barrier_width must be positive")
    barrier_left = parameters.barrier_x - parameters.barrier_width / 2.0
    barrier_right = parameters.barrier_x + parameters.barrier_width / 2.0
    if not -half_length < parameters.barrier_x < half_length:
        raise ValueError("barrier_x must lie inside the grid")
    if barrier_left <= -half_length or barrier_right >= half_length:
        raise ValueError("barrier must fit entirely inside the grid")
    if parameters.detector_buffer < 0:
        raise ValueError("detector_buffer must not be negative")
    detector_x = barrier_right + parameters.detector_buffer
    usable_right_edge = (
        half_length - parameters.boundary_width
        if parameters.absorbing_boundary
        else half_length
    )
    if detector_x >= usable_right_edge:
        raise ValueError("detector must lie before the absorbing or periodic grid edge")
    if parameters.measurement_kind == "zeno":
        if parameters.measurement_interval is None or parameters.measurement_interval <= 0:
            raise ValueError("Zeno measurements require a positive interval")
        if parameters.measurement_interval < parameters.dt:
            raise ValueError("measurement_interval must be at least dt")


def _extent(result: ExperimentResult) -> tuple[float, float, float, float]:
    grid = result.grid
    return (
        float(grid.x[0] - grid.dx / 2.0),
        float(grid.x[-1] + grid.dx / 2.0),
        float(grid.y[0] - grid.dy / 2.0),
        float(grid.y[-1] + grid.dy / 2.0),
    )
