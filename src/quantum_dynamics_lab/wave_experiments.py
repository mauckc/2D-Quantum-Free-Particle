"""Configuration adapter for scalar-optics forward experiments."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from .optics import (
    ComplexArray,
    OpticalGrid,
    apply_thin_element,
    circular_aperture,
    field_power,
    gaussian_mode,
    hermite_gaussian_mode,
    make_optical_grid,
    phase_only_mask,
    rectangular_aperture,
    thin_lens,
)
from .optics_propagation import (
    angular_spectrum_retained_power,
    propagate_angular_spectrum,
    propagate_fresnel,
)
from .wave_config import OpticalSourceConfig, OpticalStepConfig, WavePropagationConfig


@dataclass(frozen=True)
class WaveRun:
    config: WavePropagationConfig
    grid: OpticalGrid
    input_field: ComplexArray
    final_field: ComplexArray
    metrics: dict[str, Any]


def run_wave_propagation(config: WavePropagationConfig) -> WaveRun:
    """Execute ordered thin-element and propagation steps with NumPy."""

    grid = make_optical_grid(
        config.grid.shape,
        config.grid.extent_m,
        config.grid.wavelength_m,
    )
    input_field = _make_source(config.source, grid)
    field = input_field.copy()
    stage_metrics = []
    total_distance = 0.0
    for index, step in enumerate(config.steps):
        before = field_power(field, grid)
        details: dict[str, Any] = {"index": index, "kind": step.kind}
        field = _apply_step(field, grid, step, details)
        if step.kind.lower() == "propagate":
            total_distance += float(step.distance_m or 0.0)
        details["power_before"] = before
        details["power_after"] = field_power(field, grid)
        stage_metrics.append(details)

    input_power = field_power(input_field, grid)
    final_power = field_power(field, grid)
    metrics = {
        "name": config.name,
        "scope": "monochromatic coherent scalar free-space optics",
        "backend": "numpy",
        "wavelength_m": grid.wavelength,
        "shape": list(grid.shape),
        "extent_m": list(grid.extent),
        "sampling_m": [grid.dy, grid.dx],
        "input_power": input_power,
        "final_power": final_power,
        "power_efficiency": final_power / input_power,
        "total_distance_m": total_distance,
        "output_second_moment_waist_m": list(_second_moment_waist(field, grid)),
        "stages": stage_metrics,
    }
    return WaveRun(
        config=config,
        grid=grid,
        input_field=input_field,
        final_field=field,
        metrics=metrics,
    )


def _make_source(source: OpticalSourceConfig, grid: OpticalGrid) -> ComplexArray:
    kind = source.kind.lower()
    if kind == "gaussian":
        return gaussian_mode(grid, source.waist_m, center=source.center_m)
    if kind in {"hermite_gaussian", "hg"}:
        return hermite_gaussian_mode(
            grid,
            source.order,
            source.waist_m,
            center=source.center_m,
        )
    if kind in {"array", "complex_array"}:
        if source.path is None:
            raise ValueError("array source requires source.path")
        if source.path.suffix.lower() == ".npy":
            values = np.load(source.path)
        elif source.path.suffix.lower() == ".npz":
            with np.load(source.path) as archive:
                if source.array_key not in archive:
                    raise ValueError(
                        f"array source key {source.array_key!r} is not present"
                    )
                values = archive[source.array_key]
        else:
            raise ValueError("array source path must end in .npy or .npz")
        from .optics import as_complex_field

        return as_complex_field(values, grid)
    raise ValueError(f"unsupported optical source kind: {source.kind}")


def _apply_step(
    field: ComplexArray,
    grid: OpticalGrid,
    step: OpticalStepConfig,
    details: dict[str, Any],
) -> ComplexArray:
    kind = step.kind.lower()
    if kind == "lens":
        if step.focal_length_m is None:
            raise ValueError("lens step requires focal_length_m")
        details["focal_length_m"] = step.focal_length_m
        return apply_thin_element(
            field,
            thin_lens(grid, step.focal_length_m, center=step.center_m),
            grid,
        )
    if kind == "circular_aperture":
        if step.radius_m is None:
            raise ValueError("circular_aperture step requires radius_m")
        details["radius_m"] = step.radius_m
        return apply_thin_element(
            field,
            circular_aperture(grid, step.radius_m, center=step.center_m),
            grid,
        )
    if kind == "rectangular_aperture":
        if step.size_m is None:
            raise ValueError("rectangular_aperture step requires size_m")
        details["size_m"] = list(step.size_m)
        return apply_thin_element(
            field,
            rectangular_aperture(grid, step.size_m, center=step.center_m),
            grid,
        )
    if kind == "phase_mask":
        phase = _analytic_phase(step, grid)
        details.update(
            {
                "phase_profile": step.phase_profile,
                "phase_amplitude_rad": step.phase_amplitude_rad,
                "phase_period_m": step.phase_period_m,
                "axis": step.axis,
            }
        )
        return apply_thin_element(field, phase_only_mask(phase, grid), grid)
    if kind == "propagate":
        if step.distance_m is None:
            raise ValueError("propagate step requires distance_m")
        model = step.model.lower()
        details.update({"model": model, "distance_m": step.distance_m})
        if model == "fresnel":
            return propagate_fresnel(field, grid, step.distance_m)
        if model in {"angular_spectrum", "blas"}:
            details["retained_spectral_power"] = angular_spectrum_retained_power(
                field,
                grid,
                step.distance_m,
            )
            return propagate_angular_spectrum(field, grid, step.distance_m)
        raise ValueError(f"unsupported propagation model: {step.model}")
    raise ValueError(f"unsupported optical step kind: {step.kind}")


def _analytic_phase(step: OpticalStepConfig, grid: OpticalGrid) -> np.ndarray:
    if step.phase_profile.lower() != "sinusoidal":
        raise ValueError(f"unsupported phase profile: {step.phase_profile}")
    if step.phase_period_m is None or step.phase_period_m <= 0.0:
        raise ValueError("sinusoidal phase mask requires positive phase_period_m")
    if abs(step.phase_amplitude_rad) > math.pi:
        raise ValueError("phase_amplitude_rad must lie within [-pi, pi]")
    if step.axis.lower() == "x":
        coordinate = grid.X - step.center_m[1]
    elif step.axis.lower() == "y":
        coordinate = grid.Y - step.center_m[0]
    else:
        raise ValueError("phase-mask axis must be 'x' or 'y'")
    return step.phase_amplitude_rad * np.sin(
        2.0 * math.pi * coordinate / step.phase_period_m
    )


def _second_moment_waist(
    field: ComplexArray,
    grid: OpticalGrid,
) -> tuple[float, float]:
    intensity = np.abs(field) ** 2
    power = field_power(field, grid)
    mean_y = float(np.sum(intensity * grid.Y) * grid.sample_area / power)
    mean_x = float(np.sum(intensity * grid.X) * grid.sample_area / power)
    variance_y = float(
        np.sum(intensity * (grid.Y - mean_y) ** 2) * grid.sample_area / power
    )
    variance_x = float(
        np.sum(intensity * (grid.X - mean_x) ** 2) * grid.sample_area / power
    )
    return 2.0 * math.sqrt(variance_y), 2.0 * math.sqrt(variance_x)
