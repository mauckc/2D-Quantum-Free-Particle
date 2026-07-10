from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .config import ExperimentConfig
from .grid import Grid, integrate_probability, make_grid, norm, normalize
from .potentials import (
    build_potential,
    lower_slit_mask,
    transmission_mask,
    upper_slit_mask,
    which_path_masks,
)
from .solver import SplitStepSolver
from .wavepackets import gaussian_packet


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


@dataclass(frozen=True)
class ExperimentResult:
    config: ExperimentConfig
    grid: Grid
    potential: FloatArray
    masks: dict[str, BoolArray]
    arrays: dict[str, NDArray]
    metrics: dict[str, object]


def run_experiment(config: ExperimentConfig) -> ExperimentResult:
    experiment = config.experiment.lower()
    if experiment in {"free_packet", "barrier_zeno"}:
        return _run_single_wavefunction(config)
    if experiment == "double_slit":
        return _run_double_slit(config)
    raise ValueError(f"unsupported experiment: {config.experiment}")


def _run_single_wavefunction(config: ExperimentConfig) -> ExperimentResult:
    grid = make_grid(config.grid)
    potential = build_potential(config.potential, grid)
    psi = gaussian_packet(config.wave_packet, grid)
    solver = SplitStepSolver(grid, potential, config.solver)
    mask = transmission_mask(config.potential, grid)
    masks = {"transmission": mask, "non_transmission": ~mask}

    frames: list[ComplexArray] = [psi.copy()]
    times: list[float] = [0.0]
    norm_series: list[float] = [norm(psi, grid)]
    transmission_series: list[float] = [_probability_in_mask(psi, grid, mask)]
    detected_series: list[float] = [0.0]
    survival_series: list[float] = [1.0]
    survival_weight = 1.0
    detected_probability = 0.0
    next_measurement = _first_measurement_time(config)

    for step in range(1, solver.step_count + 1):
        psi = solver.step(psi)
        time = step * config.solver.dt

        if _measurement_due(config, time, next_measurement):
            transmitted = _probability_in_mask(psi, grid, mask)
            detected_probability += survival_weight * transmitted
            survival_weight *= max(0.0, 1.0 - transmitted)
            psi = np.where(mask, 0.0, psi)
            if norm(psi, grid) > 0:
                psi = normalize(psi, grid)
            next_measurement += config.measurement.interval or config.solver.tf

        norm_series.append(norm(psi, grid))
        transmission_series.append(_probability_in_mask(psi, grid, mask))
        detected_series.append(detected_probability)
        survival_series.append(survival_weight)
        if step % solver.frame_stride == 0 or step == solver.step_count:
            frames.append(psi.copy())
            times.append(time)

    arrays: dict[str, NDArray] = {
        "times": np.asarray(times, dtype=np.float64),
        "psi_frames": np.asarray(frames, dtype=np.complex128),
        "metric_times": np.arange(solver.step_count + 1, dtype=np.float64)
        * config.solver.dt,
        "norm": np.asarray(norm_series, dtype=np.float64),
        "transmission_probability": np.asarray(transmission_series, dtype=np.float64),
        "detected_probability": np.asarray(detected_series, dtype=np.float64),
        "survival_weight": np.asarray(survival_series, dtype=np.float64),
    }
    metrics = _base_metrics(config, grid, arrays)
    metrics.update(
        {
            "final_norm": float(norm_series[-1]),
            "final_transmission_probability": float(transmission_series[-1]),
            "final_detected_probability": float(detected_series[-1]),
            "final_survival_weight": float(survival_series[-1]),
        }
    )
    return ExperimentResult(config, grid, potential, masks, arrays, metrics)


def _run_double_slit(config: ExperimentConfig) -> ExperimentResult:
    grid = make_grid(config.grid)
    potential = build_potential(config.potential, grid)
    psi = gaussian_packet(config.wave_packet, grid)
    solver = SplitStepSolver(grid, potential, config.solver)
    transmission = transmission_mask(config.potential, grid)
    upper_slit = upper_slit_mask(config.potential, grid)
    lower_slit = lower_slit_mask(config.potential, grid)
    upper_path, lower_path = which_path_masks(config.potential, grid)
    masks = {
        "transmission": transmission,
        "upper_slit": upper_slit,
        "lower_slit": lower_slit,
        "upper_path": upper_path,
        "lower_path": lower_path,
    }

    coherent = psi.copy()
    upper_branch: ComplexArray | None = None
    lower_branch: ComplexArray | None = None
    measurement_time = _which_path_time(config)
    measured = False

    coherent_frames: list[ComplexArray] = [coherent.copy()]
    which_path_frames: list[FloatArray] = [np.abs(coherent) ** 2]
    times: list[float] = [0.0]
    norm_series: list[float] = [norm(coherent, grid)]
    which_path_norm_series: list[float] = [norm(coherent, grid)]

    for step in range(1, solver.step_count + 1):
        coherent = solver.step(coherent)
        time = step * config.solver.dt

        if not measured and time >= measurement_time:
            upper_branch = np.where(upper_path, coherent, 0.0).astype(np.complex128)
            lower_branch = np.where(lower_path, coherent, 0.0).astype(np.complex128)
            measured = True
        elif measured and upper_branch is not None and lower_branch is not None:
            upper_branch = solver.step(upper_branch)
            lower_branch = solver.step(lower_branch)

        if measured and upper_branch is not None and lower_branch is not None:
            which_path_probability = np.abs(upper_branch) ** 2 + np.abs(lower_branch) ** 2
        else:
            which_path_probability = np.abs(coherent) ** 2

        norm_series.append(norm(coherent, grid))
        which_path_norm_series.append(integrate_probability(which_path_probability, grid))
        if step % solver.frame_stride == 0 or step == solver.step_count:
            coherent_frames.append(coherent.copy())
            which_path_frames.append(which_path_probability.copy())
            times.append(time)

    coherent_probability = np.abs(coherent) ** 2
    contrast_coherent = _screen_contrast(coherent_probability, transmission, grid)
    contrast_which = _screen_contrast(which_path_frames[-1], transmission, grid)

    arrays: dict[str, NDArray] = {
        "times": np.asarray(times, dtype=np.float64),
        "psi_frames": np.asarray(coherent_frames, dtype=np.complex128),
        "which_path_probability_frames": np.asarray(which_path_frames, dtype=np.float64),
        "metric_times": np.arange(solver.step_count + 1, dtype=np.float64)
        * config.solver.dt,
        "norm": np.asarray(norm_series, dtype=np.float64),
        "which_path_norm": np.asarray(which_path_norm_series, dtype=np.float64),
        "coherent_screen_profile": _screen_profile(coherent_probability, transmission),
        "which_path_screen_profile": _screen_profile(which_path_frames[-1], transmission),
    }
    metrics = _base_metrics(config, grid, arrays)
    metrics.update(
        {
            "final_norm": float(norm_series[-1]),
            "final_which_path_norm": float(which_path_norm_series[-1]),
            "coherent_interference_contrast": float(contrast_coherent),
            "which_path_interference_contrast": float(contrast_which),
        }
    )
    return ExperimentResult(config, grid, potential, masks, arrays, metrics)


def _first_measurement_time(config: ExperimentConfig) -> float:
    if config.measurement.time is not None:
        return config.measurement.time
    if config.measurement.interval is None:
        return config.solver.tf + config.solver.dt
    return config.measurement.interval


def _measurement_due(config: ExperimentConfig, time: float, next_measurement: float) -> bool:
    return (
        config.measurement.kind.lower() == "zeno"
        and config.measurement.interval is not None
        and config.measurement.interval > 0
        and time + config.solver.dt / 2.0 >= next_measurement
    )


def _which_path_time(config: ExperimentConfig) -> float:
    if config.measurement.time is not None:
        return config.measurement.time
    px = config.wave_packet.momentum[0]
    if px > 0:
        travel_time = (
            config.potential.barrier_x - config.wave_packet.center[0]
        ) * config.solver.mass / px
        if 0 < travel_time < config.solver.tf:
            return travel_time
    return 0.5 * config.solver.tf


def _probability_in_mask(psi: ComplexArray, grid: Grid, mask: BoolArray) -> float:
    return integrate_probability(np.where(mask, np.abs(psi) ** 2, 0.0), grid)


def _screen_profile(probability: FloatArray, transmission: BoolArray) -> FloatArray:
    columns = np.where(np.any(transmission, axis=0))[0]
    if len(columns) == 0:
        return probability[:, -1]
    return probability[:, columns[-1]]


def _screen_contrast(probability: FloatArray, transmission: BoolArray, grid: Grid) -> float:
    profile = _screen_profile(probability, transmission)
    maximum = float(np.max(profile))
    minimum = float(np.min(profile))
    total = maximum + minimum
    if total <= 0:
        return 0.0
    return (maximum - minimum) / total


def _base_metrics(
    config: ExperimentConfig, grid: Grid, arrays: dict[str, NDArray]
) -> dict[str, object]:
    return {
        "name": config.name,
        "experiment": config.experiment,
        "grid_n": grid.n,
        "grid_length": grid.length,
        "dx": grid.dx,
        "dt": config.solver.dt,
        "tf": config.solver.tf,
        "frame_count": int(len(arrays["times"])),
    }
