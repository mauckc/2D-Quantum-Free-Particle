"""Independent model-transfer and convergence validation for the flagship."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json

import numpy as np

from .flagship import BASELINE_NAMES, FlagshipResult, run_flagship
from .flagship_config import FlagshipConfig, PerturbationScenario
from .objectives import intensity_similarity, mode_overlap, mode_power_efficiency
from .optics import circular_aperture, gaussian_mode, hermite_gaussian_mode, make_optical_grid
from .optics_propagation import (
    angular_spectrum_retained_power,
    propagate_angular_spectrum,
    propagate_fresnel,
)


DECLARED_LIMITS = {
    "maximum_overlap_transfer_loss": 1.0e-3,
    "maximum_efficiency_transfer_loss": 1.0e-3,
    "minimum_blas_retained_spectral_power": 1.0 - 1.0e-12,
    "maximum_grid_finest_pair_overlap_change": 2.0e-2,
    "maximum_padding_finest_pair_overlap_change": 5.0e-3,
    "maximum_aperture_finest_pair_overlap_change": 5.0e-3,
    "maximum_step_finest_pair_overlap_change": 1.0e-8,
    "maximum_boundary_power_fraction": 1.0e-3,
    "reproduction_absolute_tolerance": 2.0e-10,
    "reproduction_relative_tolerance": 2.0e-10,
}


@dataclass(frozen=True)
class FlagshipValidationResult:
    metrics: dict[str, Any]


def validate_flagship(
    config: FlagshipConfig,
    flagship: FlagshipResult | None = None,
    out_dir: str | Path | None = None,
) -> FlagshipValidationResult:
    """Run model-transfer and numerical convergence on fixed optimized designs."""

    if flagship is None:
        flagship = run_flagship(config)
    base_grid = make_optical_grid(
        config.grid.shape,
        config.grid.extent_m,
        config.grid.wavelength_m,
    )
    source = gaussian_mode(base_grid, config.system.waist_m)
    target = hermite_gaussian_mode(base_grid, (1, 0), config.system.waist_m)
    nominal = PerturbationScenario(name="nominal", kind="nominal")
    scenario_sequence = (nominal, *config.held_out_scenarios)
    model_transfer: dict[str, Any] = {}
    all_overlap_losses = []
    all_efficiency_losses = []
    all_retained = []
    for name in BASELINE_NAMES:
        comparisons = []
        for scenario in scenario_sequence:
            fresnel = _evaluate(
                flagship.designs[name],
                config,
                base_grid,
                source,
                target,
                scenario,
                model="fresnel",
            )
            angular = _evaluate(
                flagship.designs[name],
                config,
                base_grid,
                source,
                target,
                scenario,
                model="angular_spectrum",
            )
            overlap_loss = fresnel["mode_overlap"] - angular["mode_overlap"]
            efficiency_loss = (
                fresnel["mode_power_efficiency"]
                - angular["mode_power_efficiency"]
            )
            all_overlap_losses.append(abs(overlap_loss))
            all_efficiency_losses.append(abs(efficiency_loss))
            all_retained.append(angular["minimum_retained_spectral_power"])
            comparisons.append(
                {
                    "name": scenario.name,
                    "kind": scenario.kind,
                    "fresnel": fresnel,
                    "angular_spectrum": angular,
                    "overlap_transfer_loss": overlap_loss,
                    "efficiency_transfer_loss": efficiency_loss,
                }
            )
        model_transfer[name] = comparisons

    robust_design = flagship.designs["three_robust_masks"]
    grid_series = _grid_series(robust_design, config, base_grid)
    padding_series = _padding_series(robust_design, config, base_grid)
    aperture_series = _aperture_series(robust_design, config, base_grid)
    step_series = _step_series(robust_design, config, base_grid)
    boundary_fraction = padding_series[-1]["boundary_power_fraction"]
    convergence = {
        "grid": grid_series,
        "padding": padding_series,
        "aperture": aperture_series,
        "step": step_series,
        "grid_finest_pair_overlap_change": _finest_change(grid_series),
        "padding_finest_pair_overlap_change": _finest_change(padding_series),
        "aperture_finest_pair_overlap_change": _finest_change(aperture_series),
        "step_finest_pair_overlap_change": _finest_change(step_series),
        "boundary_power_fraction": boundary_fraction,
    }
    summary = {
        "maximum_absolute_overlap_transfer_loss": max(all_overlap_losses),
        "maximum_absolute_efficiency_transfer_loss": max(all_efficiency_losses),
        "minimum_blas_retained_spectral_power": min(all_retained),
    }
    checks = {
        "model_transfer_overlap": summary[
            "maximum_absolute_overlap_transfer_loss"
        ]
        <= DECLARED_LIMITS["maximum_overlap_transfer_loss"],
        "model_transfer_efficiency": summary[
            "maximum_absolute_efficiency_transfer_loss"
        ]
        <= DECLARED_LIMITS["maximum_efficiency_transfer_loss"],
        "blas_retained_power": summary["minimum_blas_retained_spectral_power"]
        >= DECLARED_LIMITS["minimum_blas_retained_spectral_power"],
        "grid_convergence": convergence["grid_finest_pair_overlap_change"]
        <= DECLARED_LIMITS["maximum_grid_finest_pair_overlap_change"],
        "padding_convergence": convergence["padding_finest_pair_overlap_change"]
        <= DECLARED_LIMITS["maximum_padding_finest_pair_overlap_change"],
        "aperture_convergence": convergence["aperture_finest_pair_overlap_change"]
        <= DECLARED_LIMITS["maximum_aperture_finest_pair_overlap_change"],
        "step_convergence": convergence["step_finest_pair_overlap_change"]
        <= DECLARED_LIMITS["maximum_step_finest_pair_overlap_change"],
        "boundary_power": boundary_fraction
        <= DECLARED_LIMITS["maximum_boundary_power_fraction"],
    }
    metrics = {
        "schema_version": 1,
        "name": config.name,
        "model_transfer": model_transfer,
        "model_transfer_summary": summary,
        "convergence": convergence,
        "declared_limits": DECLARED_LIMITS,
        "checks": checks,
        "trusted": all(checks.values()),
        "trusted_regime": (
            "monochromatic coherent scalar free-space fields within the committed "
            "sampling, padding, aperture, retained-spectrum, and paraxial/BLAS "
            "transfer thresholds"
        ),
        "maxwell_fdtd": {
            "performed": False,
            "required": False,
            "rationale": (
                "Optional reduced-scale full-wave validation was not performed; "
                "the millimetre-scale free-space device is outside the practical "
                "CPU CI scope and ADR 0004 makes scalar-model transfer the release gate."
            ),
        },
    }
    result = FlagshipValidationResult(metrics=metrics)
    if out_dir is not None:
        _write_validation_artifacts(Path(out_dir), result)
    return result


def _evaluate(
    phases: np.ndarray,
    config: FlagshipConfig,
    grid,
    source: np.ndarray,
    target: np.ndarray,
    scenario: PerturbationScenario,
    *,
    model: str,
    aperture_radius: float | None = None,
    step_count: int = 1,
) -> dict[str, float]:
    aperture = circular_aperture(
        grid,
        config.system.aperture_radius_m if aperture_radius is None else aperture_radius,
    )
    scenario_grid = make_optical_grid(
        grid.shape,
        grid.extent,
        grid.wavelength * scenario.wavelength_scale,
    )
    shifted = np.roll(
        phases,
        shift=scenario.alignment_pixels,
        axis=(-2, -1),
    ) * scenario.phase_depth_scale
    field = source.copy()
    retained = 1.0
    for phase, distance in zip(
        shifted,
        config.system.plane_distances_m,
        strict=True,
    ):
        field = field * aperture * np.exp(1j * phase)
        substep = distance * scenario.spacing_scale / step_count
        for _ in range(step_count):
            if model == "fresnel":
                field = propagate_fresnel(field, scenario_grid, substep)
            elif model == "angular_spectrum":
                retained = min(
                    retained,
                    angular_spectrum_retained_power(field, scenario_grid, substep),
                )
                field = propagate_angular_spectrum(field, scenario_grid, substep)
            else:
                raise ValueError(f"unsupported validation model: {model}")
    power = float(np.sum(np.abs(field) ** 2) * grid.sample_area)
    return {
        "mode_overlap": mode_overlap(field, target, grid.sample_area),
        "intensity_similarity": intensity_similarity(field, target, grid.sample_area),
        "mode_power_efficiency": mode_power_efficiency(
            field,
            target,
            1.0,
            grid.sample_area,
        ),
        "output_power": power,
        "minimum_retained_spectral_power": retained,
        "boundary_power_fraction": _boundary_power_fraction(field, grid),
    }


def _grid_series(phases: np.ndarray, config: FlagshipConfig, base_grid) -> list[dict[str, float]]:
    series = []
    for size in (24, 32, 48, 64):
        grid = make_optical_grid(size, config.grid.extent_m, config.grid.wavelength_m)
        resized = _resample_phase_masks(phases, base_grid, grid)
        source = gaussian_mode(grid, config.system.waist_m)
        target = hermite_gaussian_mode(grid, (1, 0), config.system.waist_m)
        metrics = _evaluate(
            resized,
            config,
            grid,
            source,
            target,
            PerturbationScenario("nominal", "nominal"),
            model="fresnel",
        )
        series.append({"grid_size": float(size), **metrics})
    return series


def _padding_series(
    phases: np.ndarray,
    config: FlagshipConfig,
    base_grid,
) -> list[dict[str, float]]:
    series = []
    for factor in (1.0, 1.5, 2.0):
        shape = tuple(int(round(value * factor)) for value in base_grid.shape)
        extent = tuple(value * factor for value in base_grid.extent)
        grid = make_optical_grid(shape, extent, config.grid.wavelength_m)
        resized = _resample_phase_masks(phases, base_grid, grid)
        source = gaussian_mode(grid, config.system.waist_m)
        target = hermite_gaussian_mode(grid, (1, 0), config.system.waist_m)
        metrics = _evaluate(
            resized,
            config,
            grid,
            source,
            target,
            PerturbationScenario("nominal", "nominal"),
            model="fresnel",
        )
        series.append({"padding_factor": factor, **metrics})
    return series


def _aperture_series(
    phases: np.ndarray,
    config: FlagshipConfig,
    grid,
) -> list[dict[str, float]]:
    source = gaussian_mode(grid, config.system.waist_m)
    target = hermite_gaussian_mode(grid, (1, 0), config.system.waist_m)
    series = []
    for radius in (1.2e-3, 1.5e-3, 1.8e-3, 1.95e-3):
        metrics = _evaluate(
            phases,
            config,
            grid,
            source,
            target,
            PerturbationScenario("nominal", "nominal"),
            model="fresnel",
            aperture_radius=radius,
        )
        series.append({"aperture_radius_m": radius, **metrics})
    return series


def _step_series(phases: np.ndarray, config: FlagshipConfig, grid) -> list[dict[str, float]]:
    source = gaussian_mode(grid, config.system.waist_m)
    target = hermite_gaussian_mode(grid, (1, 0), config.system.waist_m)
    series = []
    for step_count in (1, 2, 4):
        metrics = _evaluate(
            phases,
            config,
            grid,
            source,
            target,
            PerturbationScenario("nominal", "nominal"),
            model="fresnel",
            step_count=step_count,
        )
        series.append({"substeps_per_plane": float(step_count), **metrics})
    return series


def _resample_phase_masks(phases: np.ndarray, old_grid, new_grid) -> np.ndarray:
    if old_grid.shape == new_grid.shape and old_grid.extent == new_grid.extent:
        return phases.copy()
    result = []
    for phase in phases:
        transmission = np.exp(1j * phase)
        intermediate = np.empty((old_grid.shape[0], new_grid.shape[1]), np.complex128)
        for row in range(old_grid.shape[0]):
            intermediate[row] = np.interp(
                new_grid.x,
                old_grid.x,
                transmission[row].real,
                left=1.0,
                right=1.0,
            ) + 1j * np.interp(
                new_grid.x,
                old_grid.x,
                transmission[row].imag,
                left=0.0,
                right=0.0,
            )
        resized = np.empty(new_grid.shape, np.complex128)
        for column in range(new_grid.shape[1]):
            resized[:, column] = np.interp(
                new_grid.y,
                old_grid.y,
                intermediate[:, column].real,
                left=1.0,
                right=1.0,
            ) + 1j * np.interp(
                new_grid.y,
                old_grid.y,
                intermediate[:, column].imag,
                left=0.0,
                right=0.0,
            )
        result.append(np.angle(resized))
    return np.asarray(result, dtype=np.float64)


def _boundary_power_fraction(field: np.ndarray, grid, width: int = 2) -> float:
    intensity = np.abs(field) ** 2
    mask = np.zeros(grid.shape, dtype=bool)
    mask[:width, :] = True
    mask[-width:, :] = True
    mask[:, :width] = True
    mask[:, -width:] = True
    total = float(np.sum(intensity))
    return 0.0 if total == 0.0 else float(np.sum(intensity[mask]) / total)


def _finest_change(series: list[dict[str, float]]) -> float:
    return abs(series[-1]["mode_overlap"] - series[-2]["mode_overlap"])


def _write_validation_artifacts(
    out_dir: Path,
    result: FlagshipValidationResult,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "validation_metrics.json").write_text(
        json.dumps(result.metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "model_transfer_report.md").write_text(
        _model_transfer_report(result.metrics),
        encoding="utf-8",
    )
    flagship_report = out_dir / "report.md"
    if flagship_report.exists():
        original = flagship_report.read_text(encoding="utf-8")
        flagship_report.write_text(
            original
            + "\n## Independent validation\n\n"
            + f"Trusted: **{result.metrics['trusted']}**. See "
            + "`model_transfer_report.md` and `validation_metrics.json` for "
            + "the full BLAS and convergence evidence.\n",
            encoding="utf-8",
        )


def _model_transfer_report(metrics: dict[str, Any]) -> str:
    summary = metrics["model_transfer_summary"]
    convergence = metrics["convergence"]
    checks = "\n".join(
        f"- {name}: {'PASS' if passed else 'FAIL'}"
        for name, passed in metrics["checks"].items()
    )
    return f"""# Flagship model-transfer and convergence report

Trusted: **{metrics['trusted']}**

- Maximum absolute Fresnel/BLAS overlap transfer loss: {summary['maximum_absolute_overlap_transfer_loss']:.8g}
- Maximum absolute Fresnel/BLAS efficiency transfer loss: {summary['maximum_absolute_efficiency_transfer_loss']:.8g}
- Minimum BLAS retained spectral power: {summary['minimum_blas_retained_spectral_power']:.12f}
- Grid finest-pair overlap change: {convergence['grid_finest_pair_overlap_change']:.8g}
- Padding finest-pair overlap change: {convergence['padding_finest_pair_overlap_change']:.8g}
- Aperture finest-pair overlap change: {convergence['aperture_finest_pair_overlap_change']:.8g}
- Propagation-substep finest-pair overlap change: {convergence['step_finest_pair_overlap_change']:.8g}
- Padded-window boundary power fraction: {convergence['boundary_power_fraction']:.8g}

## Declared checks

{checks}

## Trusted regime

{metrics['trusted_regime']}.

## Full-wave disposition

Maxwell/FDTD performed: **False**. {metrics['maxwell_fdtd']['rationale']}
"""
