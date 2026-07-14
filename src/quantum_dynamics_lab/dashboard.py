from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
import json
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

from .artifacts import save_run, write_metrics
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
from .render import render_comparison, render_run
from .validation import load_validation_config, run_validation
from .wave_artifacts import save_wave_run
from .wave_config import load_wave_config
from .wave_experiments import WaveRun, run_wave_propagation


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
    slit_width: float = 0.45
    slit_separation: float = 1.2
    measurement_kind: str = "zeno"
    measurement_interval: float | None = 0.12
    measurement_time: float | None = None
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


@dataclass(frozen=True)
class DashboardSweep:
    out_dir: Path
    run_paths: tuple[Path, ...]
    rows: tuple[dict[str, float], ...]
    outputs: dict[str, Path]


@dataclass(frozen=True)
class WaveDashboardRun:
    run: WaveRun
    run_path: Path


@dataclass(frozen=True)
class FlagshipDashboardSummary:
    comparison: dict[str, Any]
    validation: dict[str, Any]


def run_dashboard_wave_forward(
    config_path: str | Path,
    out_dir: str | Path,
) -> WaveDashboardRun:
    """Run a NumPy forward-optics config without importing optional JAX."""

    run = run_wave_propagation(load_wave_config(config_path))
    return WaveDashboardRun(run=run, run_path=save_wave_run(run, out_dir))


def plot_wave_field_panels(run: WaveRun) -> tuple[Figure, Figure]:
    """Plot input/output intensity and phase for dashboard display."""

    figures = []
    extent = [
        float(run.grid.x[0]),
        float(run.grid.x[-1]),
        float(run.grid.y[0]),
        float(run.grid.y[-1]),
    ]
    for label, field in (
        ("Input field", run.input_field),
        ("Output field", run.final_field),
    ):
        figure, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)
        intensity = axes[0].imshow(
            np.abs(field) ** 2,
            origin="lower",
            extent=extent,
            cmap="magma",
        )
        axes[0].set_title(f"{label} intensity")
        figure.colorbar(intensity, ax=axes[0], shrink=0.8)
        phase = axes[1].imshow(
            np.angle(field),
            origin="lower",
            extent=extent,
            cmap="twilight",
            vmin=-np.pi,
            vmax=np.pi,
        )
        axes[1].set_title(f"{label} phase")
        figure.colorbar(phase, ax=axes[1], shrink=0.8)
        for axis in axes:
            axis.set_xlabel("x")
            axis.set_ylabel("y")
        figures.append(figure)
    return tuple(figures)


def load_flagship_dashboard_summary(
    repository_root: str | Path | None = None,
) -> FlagshipDashboardSummary:
    """Load committed M4 metrics for read-only UI use without JAX."""

    root = (
        Path(repository_root)
        if repository_root is not None
        else Path(__file__).resolve().parents[2]
    )
    reference = root / "benchmarks" / "reference"
    comparison = json.loads(
        (reference / "flagship_m4.json").read_text(encoding="utf-8")
    )
    validation = json.loads(
        (reference / "flagship_validation_m4.json").read_text(encoding="utf-8")
    )
    return FlagshipDashboardSummary(comparison=comparison, validation=validation)


def build_experiment_config(parameters: DashboardParameters) -> ExperimentConfig:
    _validate_parameters(parameters)
    is_barrier = parameters.experiment == "barrier_zeno"
    is_double_slit = parameters.experiment == "double_slit"
    potential_kind = (
        "double_slit" if is_double_slit else "barrier" if is_barrier else "free"
    )
    measurement_kind = (
        "which_path" if is_double_slit else parameters.measurement_kind if is_barrier else "none"
    )
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
            kind=potential_kind,
            height=parameters.barrier_height if is_barrier or is_double_slit else 0.0,
            barrier_x=parameters.barrier_x,
            barrier_width=parameters.barrier_width,
            slit_width=parameters.slit_width,
            slit_separation=parameters.slit_separation,
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
            time=parameters.measurement_time if is_double_slit else None,
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


def run_dashboard_zeno_sweep(
    parameters: DashboardParameters,
    heights: tuple[float, ...],
    intervals: tuple[float, ...],
    out_dir: str | Path,
) -> DashboardSweep:
    if not heights or not intervals:
        raise ValueError("Zeno sweeps require at least one height and interval")
    if len(heights) * len(intervals) > 25:
        raise ValueError("interactive Zeno sweeps are limited to 25 variants")
    if any(height < 0 for height in heights):
        raise ValueError("barrier heights must not be negative")
    if any(interval < 0 for interval in intervals):
        raise ValueError("measurement intervals must not be negative")
    if not any(interval == 0 for interval in intervals):
        raise ValueError("Zeno sweeps must include an unmeasured interval of 0")
    if not any(interval > 0 for interval in intervals):
        raise ValueError("Zeno sweeps require at least one measured interval")

    out_dir = Path(out_dir)
    run_paths: list[Path] = []
    run_dirs: list[Path] = []
    rows: list[dict[str, float]] = []
    variant_metrics: list[dict[str, object]] = []
    for height in sorted(set(heights)):
        for interval in sorted(set(intervals)):
            measured = interval > 0
            name = _sweep_variant_name(height, interval)
            variant = replace(
                parameters,
                name=name,
                experiment="barrier_zeno",
                barrier_height=height,
                measurement_kind="zeno" if measured else "none",
                measurement_interval=interval if measured else None,
                measurement_time=None,
            )
            result = run_experiment(build_experiment_config(variant))
            run_dir = out_dir / "runs" / name
            run_path = save_run(result, run_dir)
            run_paths.append(run_path)
            run_dirs.append(run_dir)
            variant_metrics.append(result.metrics)
            rows.append(
                {
                    "barrier_height": float(height),
                    "measurement_interval": float(interval),
                    "transmission": float(
                        result.metrics["final_unconditional_transmission_probability"]
                    ),
                    "detected": float(result.metrics["final_detected_probability"]),
                    "survival": float(result.metrics["final_survival_weight"]),
                    "absorbed": float(result.metrics["final_absorbed_probability"]),
                }
            )

    render_comparison(run_dirs, out_dir)
    write_metrics(
        out_dir / "metrics.json",
        {"name": "dashboard_zeno_sweep", "variants": variant_metrics},
    )
    output_names = {
        "comparison": "comparison.png",
        "heatmap": "zeno_transmission_heatmap.png",
        "csv": "comparison.csv",
        "report": "report.md",
        "html_report": "report.html",
        "explorer": "explorer.html",
    }
    outputs = {
        key: out_dir / filename
        for key, filename in output_names.items()
        if (out_dir / filename).exists()
    }
    return DashboardSweep(
        out_dir=out_dir,
        run_paths=tuple(run_paths),
        rows=tuple(rows),
        outputs=outputs,
    )


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
    summary = {
        "final_norm": float(result.metrics["final_norm"]),
        "max_norm_drift": float(np.max(np.abs(result.arrays["norm"] - 1.0))),
    }
    if result.config.experiment == "double_slit":
        summary.update(
            {
                "final_which_path_norm": float(result.metrics["final_which_path_norm"]),
                "coherent_interference_contrast": float(
                    result.metrics["coherent_interference_contrast"]
                ),
                "which_path_interference_contrast": float(
                    result.metrics["which_path_interference_contrast"]
                ),
            }
        )
    else:
        summary.update(
            {
                "final_unconditional_transmission_probability": float(
                    result.metrics["final_unconditional_transmission_probability"]
                ),
                "final_detected_probability": float(
                    result.metrics["final_detected_probability"]
                ),
                "final_survival_weight": float(result.metrics["final_survival_weight"]),
                "final_absorbed_probability": float(
                    result.metrics["final_absorbed_probability"]
                ),
            }
        )
    return summary


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
    probability_axis.set_title(
        "Coherent probability density"
        if result.config.experiment == "double_slit"
        else "Probability density"
    )
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


def plot_double_slit_panels(
    result: ExperimentResult, frame_index: int
) -> tuple[Figure, ...]:
    if result.config.experiment != "double_slit":
        raise ValueError("double-slit panels require a double_slit result")
    which_path_frames = result.arrays["which_path_probability_frames"]
    if not 0 <= frame_index < len(which_path_frames):
        raise IndexError("frame_index is outside the saved frame range")

    density_figure, density_axis = plt.subplots(figsize=(5.4, 4.0))
    density_image = density_axis.imshow(
        which_path_frames[frame_index],
        origin="lower",
        extent=_extent(result),
        cmap="cividis",
        aspect="equal",
    )
    density_axis.set_title("Which-path probability density")
    density_axis.set_xlabel("x")
    density_axis.set_ylabel("y")
    density_figure.colorbar(density_image, ax=density_axis, label="probability density")

    profile_figure, profile_axis = plt.subplots(figsize=(5.4, 4.0))
    profile_axis.plot(
        result.grid.y,
        result.arrays["coherent_screen_profile"],
        label="coherent",
    )
    profile_axis.plot(
        result.grid.y,
        result.arrays["which_path_screen_profile"],
        label="which-path",
    )
    profile_axis.set_title("Final screen profile")
    profile_axis.set_xlabel("y")
    profile_axis.set_ylabel("integrated probability")
    profile_axis.grid(True, alpha=0.25)
    profile_axis.legend(loc="best")
    density_figure.tight_layout()
    profile_figure.tight_layout()
    return density_figure, profile_figure


def plot_diagnostic_panels(result: ExperimentResult) -> tuple[Figure, ...]:
    arrays = result.arrays
    metric_times = arrays["metric_times"]
    norm_figure, norm_axis = plt.subplots(figsize=(5.4, 3.6))
    probability_figure, probability_axis = plt.subplots(figsize=(5.4, 3.6))

    norm_drift = arrays["norm"] - 1.0
    if result.config.experiment == "double_slit":
        norm_axis.plot(
            metric_times, norm_drift, color="#2563eb", label="coherent norm - 1"
        )
        norm_axis.plot(
            metric_times,
            arrays["which_path_norm"] - 1.0,
            color="#b45309",
            label="which-path norm - 1",
        )
        norm_axis.set_title("Branch norm change with absorbing edges")
    else:
        norm_axis.plot(
            metric_times, norm_drift, color="#2563eb", label="conditional norm - 1"
        )
        norm_axis.set_title("Norm drift")
    norm_axis.axhline(0.0, color="#6b7280", linewidth=0.8)
    norm_axis.set_xlabel("time")
    norm_axis.set_ylabel("probability mass")
    norm_axis.grid(True, alpha=0.25)
    norm_axis.legend(loc="best")

    if result.config.experiment == "double_slit":
        labels = ("coherent", "which-path")
        values = (
            float(result.metrics["coherent_interference_contrast"]),
            float(result.metrics["which_path_interference_contrast"]),
        )
        probability_axis.bar(labels, values, color=("#2563eb", "#b45309"))
        probability_axis.set_title("Interference contrast")
        probability_axis.set_ylabel("contrast")
        probability_axis.set_ylim(0.0, max(1.0, max(values) * 1.15))
        probability_axis.grid(True, axis="y", alpha=0.25)
        for index, value in enumerate(values):
            probability_axis.text(index, value, f"{value:.3f}", ha="center", va="bottom")
    elif result.config.experiment == "barrier_zeno":
        series = (
            ("unconditional_transmission_probability", "no-click transmission"),
            ("detected_probability", "detector clicks"),
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
    if parameters.experiment not in {"free_packet", "barrier_zeno", "double_slit"}:
        raise ValueError("unsupported dashboard experiment")
    if parameters.measurement_kind not in {"none", "zeno", "which_path"}:
        raise ValueError("unsupported measurement kind")
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
    if parameters.experiment == "free_packet":
        return
    if parameters.barrier_height < 0:
        raise ValueError("barrier_height must not be negative")
    if parameters.barrier_width <= 0:
        raise ValueError("barrier_width must be positive")
    barrier_left = parameters.barrier_x - parameters.barrier_width / 2.0
    barrier_right = parameters.barrier_x + parameters.barrier_width / 2.0
    if not -half_length < parameters.barrier_x < half_length:
        raise ValueError("barrier_x must lie inside the grid")
    if barrier_left <= -half_length or barrier_right >= half_length:
        raise ValueError("barrier must fit entirely inside the grid")
    if parameters.experiment == "double_slit":
        if parameters.slit_width <= 0 or parameters.slit_separation <= 0:
            raise ValueError("slit width and separation must be positive")
        if parameters.slit_separation <= parameters.slit_width:
            raise ValueError("slit separation must exceed slit width")
        slit_outer_edge = (parameters.slit_separation + parameters.slit_width) / 2.0
        if slit_outer_edge >= half_length:
            raise ValueError("both slits must fit entirely inside the grid")
        if parameters.measurement_time is not None and not (
            0 < parameters.measurement_time < parameters.tf
        ):
            raise ValueError("which-path measurement time must lie between 0 and tf")
        return
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
        if parameters.measurement_interval > parameters.tf:
            raise ValueError("measurement_interval must not exceed tf")


def _sweep_variant_name(height: float, interval: float) -> str:
    height_label = f"{height:g}".replace("-", "m").replace(".", "p")
    if interval == 0:
        return f"h{height_label}_unmeasured"
    interval_label = f"{interval:g}".replace("-", "m").replace(".", "p")
    return f"h{height_label}_dt{interval_label}"


def _extent(result: ExperimentResult) -> tuple[float, float, float, float]:
    grid = result.grid
    return (
        float(grid.x[0] - grid.dx / 2.0),
        float(grid.x[-1] + grid.dx / 2.0),
        float(grid.y[0] - grid.dy / 2.0),
        float(grid.y[-1] + grid.dy / 2.0),
    )
