"""Robust multi-plane Gaussian-to-HG10 flagship comparison."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any
import json
import platform
import subprocess

import jax
import jax.numpy as jnp
import numpy as np

from .constraints import phase_total_variation
from .flagship_config import FlagshipConfig, PerturbationScenario
from .jax_constraints import parameterize_phase_jax, phase_total_variation_jax
from .jax_objectives import (
    intensity_similarity_jax,
    mode_overlap_jax,
    mode_power_efficiency_jax,
)
from .jax_propagation import fresnel_transfer_jax, propagate_fresnel_jax
from .optics import circular_aperture, gaussian_mode, hermite_gaussian_mode, make_optical_grid


BASELINE_NAMES = (
    "no_masks",
    "one_optimized_mask",
    "three_nominal_masks",
    "three_robust_masks",
)


@dataclass(frozen=True)
class FlagshipResult:
    metrics: dict[str, Any]
    designs: dict[str, np.ndarray]
    histories: dict[str, list[dict[str, float]]]


def run_flagship(config: FlagshipConfig, out_dir: str | Path | None = None) -> FlagshipResult:
    """Optimize the required baseline ladder and evaluate held-out scenarios."""

    jax.config.update("jax_enable_x64", True)
    _validate_config(config)
    grid = make_optical_grid(
        config.grid.shape,
        config.grid.extent_m,
        config.grid.wavelength_m,
    )
    source = gaussian_mode(grid, config.system.waist_m)
    target = hermite_gaussian_mode(grid, (1, 0), config.system.waist_m)
    aperture = circular_aperture(grid, config.system.aperture_radius_m)
    nominal = PerturbationScenario(name="nominal", kind="nominal")

    one_phase, one_history = _optimize_design(
        config,
        grid,
        source,
        target,
        aperture,
        mask_count=1,
        scenarios=(nominal,),
        seed=config.seed,
    )
    nominal_phase, nominal_history = _optimize_design(
        config,
        grid,
        source,
        target,
        aperture,
        mask_count=3,
        scenarios=(nominal,),
        seed=config.seed,
    )
    robust_phase, robust_history = _optimize_design(
        config,
        grid,
        source,
        target,
        aperture,
        mask_count=3,
        scenarios=config.optimization_scenarios,
        seed=config.seed,
    )
    zero_phase = np.zeros((3, *grid.shape), dtype=np.float64)
    one_expanded = np.concatenate(
        [one_phase, np.zeros((2, *grid.shape), dtype=np.float64)],
        axis=0,
    )
    designs = {
        "no_masks": zero_phase,
        "one_optimized_mask": one_expanded,
        "three_nominal_masks": nominal_phase,
        "three_robust_masks": robust_phase,
    }
    baselines = {
        name: _evaluate_design(
            phase_masks,
            config,
            grid,
            source,
            target,
            aperture,
            nominal,
            config.held_out_scenarios,
        )
        for name, phase_masks in designs.items()
    }
    robust_advantage = (
        baselines["three_robust_masks"]["held_out_worst"]["mode_overlap"]
        - baselines["three_nominal_masks"]["held_out_worst"]["mode_overlap"]
    )
    metrics = {
        "schema_version": 1,
        "name": config.name,
        "seed": config.seed,
        "scope": "monochromatic coherent scalar free-space Fresnel flagship",
        "required_baselines": list(BASELINE_NAMES),
        "common_physics": {
            "shape": list(grid.shape),
            "extent_m": list(grid.extent),
            "wavelength_m": grid.wavelength,
            "waist_m": config.system.waist_m,
            "plane_distances_m": list(config.system.plane_distances_m),
            "total_length_m": sum(config.system.plane_distances_m),
            "aperture_radius_m": config.system.aperture_radius_m,
            "source": "Gaussian",
            "target": "HG10",
        },
        "optimization_scenarios": [
            _scenario_dict(scenario) for scenario in config.optimization_scenarios
        ],
        "held_out_scenarios": [
            _scenario_dict(scenario) for scenario in config.held_out_scenarios
        ],
        "scenario_sets_disjoint": _scenario_sets_disjoint(config),
        "baselines": baselines,
        "robust_worst_case_overlap_advantage": robust_advantage,
        "robustness_claim_passed": robust_advantage > 0.0,
        "declared_limits": {
            "minimum_robust_worst_case_advantage": 0.01,
            "reproduction_absolute_tolerance": 2.0e-10,
            "reproduction_relative_tolerance": 2.0e-10,
        },
    }
    histories = {
        "one_optimized_mask": one_history,
        "three_nominal_masks": nominal_history,
        "three_robust_masks": robust_history,
    }
    result = FlagshipResult(metrics=metrics, designs=designs, histories=histories)
    if out_dir is not None:
        _write_flagship_artifacts(config, Path(out_dir), result)
    return result


def _optimize_design(
    config: FlagshipConfig,
    grid,
    source: np.ndarray,
    target: np.ndarray,
    aperture: np.ndarray,
    *,
    mask_count: int,
    scenarios: tuple[PerturbationScenario, ...],
    seed: int,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    source_jax = jnp.asarray(source)
    target_jax = jnp.asarray(target)
    aperture_jax = jnp.asarray(aperture)
    runtime = [
        (
            scenario,
            jnp.stack(
                [
                    fresnel_transfer_jax(
                        jnp.asarray(grid.FX),
                        jnp.asarray(grid.FY),
                        grid.wavelength * scenario.wavelength_scale,
                        distance * scenario.spacing_scale,
                    )
                    for distance in config.system.plane_distances_m
                ]
            ),
        )
        for scenario in scenarios
    ]

    def constrained(parameters):
        phases = parameterize_phase_jax(
            parameters,
            bounds=config.constraints.phase_bounds_rad,
            smoothing_passes=config.constraints.smoothing_passes,
        )
        if mask_count == 1:
            phases = jnp.concatenate(
                [phases, jnp.zeros((2, *grid.shape), dtype=phases.dtype)],
                axis=0,
            )
        return phases

    def scenario_metrics(phases, scenario, transfers):
        shifted = jnp.roll(
            phases,
            shift=scenario.alignment_pixels,
            axis=(-2, -1),
        )
        scaled = shifted * scenario.phase_depth_scale
        field = source_jax
        for plane in range(3):
            field = field * aperture_jax * jnp.exp(1j * scaled[plane])
            field = propagate_fresnel_jax(field, transfers[plane])
        return jnp.asarray(
            [
                mode_overlap_jax(field, target_jax, grid.sample_area),
                intensity_similarity_jax(field, target_jax, grid.sample_area),
                mode_power_efficiency_jax(field, target_jax, 1.0, grid.sample_area),
            ]
        )

    def objective(parameters):
        phases = constrained(parameters)
        all_metrics = jnp.stack(
            [
                scenario_metrics(phases, scenario, transfers)
                for scenario, transfers in runtime
            ]
        )
        scenario_losses = (1.0 - all_metrics[:, 0]) + config.objective.intensity_weight * (
            1.0 - all_metrics[:, 1]
        )
        aggregate = jnp.mean(scenario_losses)
        if len(runtime) > 1:
            aggregate = aggregate + config.objective.worst_case_weight * jnp.max(
                scenario_losses
            )
        active_phases = phases[:mask_count]
        tv = jnp.mean(jax.vmap(phase_total_variation_jax)(active_phases))
        loss = aggregate + config.constraints.total_variation_weight * tv
        auxiliary = jnp.asarray(
            [
                jnp.mean(all_metrics[:, 0]),
                jnp.min(all_metrics[:, 0]),
                jnp.mean(all_metrics[:, 2]),
                jnp.min(all_metrics[:, 2]),
                tv,
            ]
        )
        return loss, auxiliary

    @jax.jit
    def step(parameters, first_moment, second_moment, iteration):
        (_, _), gradient = jax.value_and_grad(objective, has_aux=True)(parameters)
        first_moment = (
            config.optimizer.beta1 * first_moment
            + (1.0 - config.optimizer.beta1) * gradient
        )
        second_moment = (
            config.optimizer.beta2 * second_moment
            + (1.0 - config.optimizer.beta2) * gradient**2
        )
        first_unbiased = first_moment / (1.0 - config.optimizer.beta1**iteration)
        second_unbiased = second_moment / (1.0 - config.optimizer.beta2**iteration)
        parameters = parameters - config.optimizer.learning_rate * first_unbiased / (
            jnp.sqrt(second_unbiased) + config.optimizer.epsilon
        )
        loss, auxiliary = objective(parameters)
        return parameters, first_moment, second_moment, loss, auxiliary

    rng = np.random.default_rng(seed)
    parameters = jnp.asarray(
        config.optimizer.initial_scale
        * rng.standard_normal((mask_count, *grid.shape)),
        dtype=jnp.float64,
    )
    first_moment = jnp.zeros_like(parameters)
    second_moment = jnp.zeros_like(parameters)
    initial_loss, initial_auxiliary = jax.jit(objective)(parameters)
    history = [_flagship_history(0, initial_loss, initial_auxiliary)]
    best_loss = float(initial_loss)
    best_parameters = parameters
    for iteration in range(1, config.optimizer.iterations + 1):
        parameters, first_moment, second_moment, loss, auxiliary = step(
            parameters,
            first_moment,
            second_moment,
            jnp.asarray(iteration, dtype=jnp.float64),
        )
        history.append(_flagship_history(iteration, loss, auxiliary))
        if float(loss) < best_loss:
            best_loss = float(loss)
            best_parameters = parameters
    best_phase = np.asarray(constrained(best_parameters))[:mask_count]
    return best_phase, history


def _flagship_history(iteration: int, loss: Any, auxiliary: Any) -> dict[str, float]:
    mean_overlap, worst_overlap, mean_efficiency, worst_efficiency, tv = np.asarray(
        auxiliary
    )
    return {
        "iteration": float(iteration),
        "loss": float(loss),
        "mean_overlap": float(mean_overlap),
        "worst_overlap": float(worst_overlap),
        "mean_efficiency": float(mean_efficiency),
        "worst_efficiency": float(worst_efficiency),
        "total_variation": float(tv),
    }


def _evaluate_design(
    phase_masks: np.ndarray,
    config: FlagshipConfig,
    grid,
    source: np.ndarray,
    target: np.ndarray,
    aperture: np.ndarray,
    nominal: PerturbationScenario,
    held_out: tuple[PerturbationScenario, ...],
) -> dict[str, Any]:
    scenarios = (nominal, *held_out)
    evaluated = [
        _evaluate_scenario(
            phase_masks,
            config,
            grid,
            source,
            target,
            aperture,
            scenario,
        )
        for scenario in scenarios
    ]
    held = evaluated[1:]
    return {
        "nominal": evaluated[0],
        "held_out": held,
        "held_out_mean": {
            key: float(np.mean([item[key] for item in held]))
            for key in ("mode_overlap", "intensity_similarity", "mode_power_efficiency")
        },
        "held_out_worst": {
            key: float(np.min([item[key] for item in held]))
            for key in ("mode_overlap", "intensity_similarity", "mode_power_efficiency")
        },
        "phase_total_variation": float(
            np.mean([phase_total_variation(phase) for phase in phase_masks])
        ),
    }


def _evaluate_scenario(
    phase_masks: np.ndarray,
    config: FlagshipConfig,
    grid,
    source: np.ndarray,
    target: np.ndarray,
    aperture: np.ndarray,
    scenario: PerturbationScenario,
) -> dict[str, Any]:
    from .objectives import intensity_similarity, mode_overlap, mode_power_efficiency
    from .optics_propagation import propagate_fresnel

    phases = np.roll(
        phase_masks,
        shift=scenario.alignment_pixels,
        axis=(-2, -1),
    ) * scenario.phase_depth_scale
    scenario_grid = make_optical_grid(
        grid.shape,
        grid.extent,
        grid.wavelength * scenario.wavelength_scale,
    )
    field = source.copy()
    for phase, distance in zip(
        phases,
        config.system.plane_distances_m,
        strict=True,
    ):
        field = field * aperture * np.exp(1j * phase)
        field = propagate_fresnel(
            field,
            scenario_grid,
            distance * scenario.spacing_scale,
        )
    return {
        "name": scenario.name,
        "kind": scenario.kind,
        "mode_overlap": mode_overlap(field, target, grid.sample_area),
        "intensity_similarity": intensity_similarity(field, target, grid.sample_area),
        "mode_power_efficiency": mode_power_efficiency(
            field,
            target,
            1.0,
            grid.sample_area,
        ),
    }


def _scenario_dict(scenario: PerturbationScenario) -> dict[str, Any]:
    return {
        "name": scenario.name,
        "kind": scenario.kind,
        "wavelength_scale": scenario.wavelength_scale,
        "alignment_pixels": list(scenario.alignment_pixels),
        "phase_depth_scale": scenario.phase_depth_scale,
        "spacing_scale": scenario.spacing_scale,
    }


def _scenario_sets_disjoint(config: FlagshipConfig) -> bool:
    optimization = {scenario.signature for scenario in config.optimization_scenarios}
    held_out = {scenario.signature for scenario in config.held_out_scenarios}
    return optimization.isdisjoint(held_out)


def _validate_config(config: FlagshipConfig) -> None:
    if tuple(sorted(set(BASELINE_NAMES))) != tuple(sorted(BASELINE_NAMES)):
        raise AssertionError("baseline names must be unique")
    if not config.optimization_scenarios or not config.held_out_scenarios:
        raise ValueError("optimization and held-out scenario sets must be nonempty")
    if not _scenario_sets_disjoint(config):
        raise ValueError("optimization and held-out scenarios must be disjoint")
    required_kinds = {"wavelength", "alignment", "phase_depth", "plane_spacing"}
    for label, scenarios in (
        ("optimization", config.optimization_scenarios),
        ("held-out", config.held_out_scenarios),
    ):
        kinds = {scenario.kind for scenario in scenarios}
        missing = required_kinds - kinds
        if missing:
            raise ValueError(f"{label} scenarios missing perturbation kinds: {missing}")
    names = [
        scenario.name
        for scenario in (*config.optimization_scenarios, *config.held_out_scenarios)
    ]
    if len(names) != len(set(names)):
        raise ValueError("scenario names must be unique")
    if config.optimizer.iterations < 1 or config.optimizer.learning_rate <= 0.0:
        raise ValueError("flagship optimizer settings must be positive")


def _write_flagship_artifacts(
    config: FlagshipConfig,
    out_dir: Path,
    result: FlagshipResult,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_dir / "designs.npz", **result.designs)
    for name, phase_masks in result.designs.items():
        np.savez_compressed(out_dir / f"{name}.npz", phase_masks=phase_masks)
    (out_dir / "metrics.json").write_text(
        json.dumps(result.metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "histories.json").write_text(
        json.dumps(result.histories, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "resolved_config.json").write_text(
        json.dumps(config.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    scenarios = {
        "optimization": result.metrics["optimization_scenarios"],
        "held_out": result.metrics["held_out_scenarios"],
    }
    (out_dir / "scenarios.json").write_text(
        json.dumps(scenarios, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    environment = {
        "created_at": datetime.now(UTC).isoformat(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "jax_version": jax.__version__,
        "dependencies": {
            name: _package_version(name)
            for name in ("jaxlib", "matplotlib", "numpy", "pytest", "scipy")
        },
        "jax_x64": bool(jax.config.x64_enabled),
        "device": str(jax.devices()[0]),
        "git_commit": _git_value("rev-parse", "HEAD"),
        "git_dirty": bool(_git_value("status", "--porcelain")),
        "seed": config.seed,
    }
    (out_dir / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "report.md").write_text(
        _flagship_report(result.metrics),
        encoding="utf-8",
    )


def _flagship_report(metrics: dict[str, Any]) -> str:
    rows = []
    for name in BASELINE_NAMES:
        baseline = metrics["baselines"][name]
        rows.append(
            f"| {name} | {baseline['nominal']['mode_overlap']:.6f} | "
            f"{baseline['held_out_mean']['mode_overlap']:.6f} | "
            f"{baseline['held_out_worst']['mode_overlap']:.6f} | "
            f"{baseline['held_out_worst']['mode_power_efficiency']:.6f} |"
        )
    claim = "passed" if metrics["robustness_claim_passed"] else "did not pass"
    return f"""# Robust Gaussian-to-HG10 flagship

| Baseline | Nominal overlap | Held-out mean overlap | Held-out worst overlap | Held-out worst efficiency |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(rows)}

Optimization and held-out scenarios are explicitly disjoint and both cover
wavelength, alignment, phase-depth, and plane-spacing perturbations. The robust
worst-case advantage over the nominal three-mask design is
`{metrics['robust_worst_case_overlap_advantage']:.6f}`; the declared robustness
comparison **{claim}** under the Fresnel model.

This is an initial scalar-Fresnel comparison. Independent BLAS transfer and
grid/padding/aperture/step convergence are separate M4 gates and must pass before
this result is presented as the v1.0 flagship conclusion.
"""


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _git_value(*args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()
