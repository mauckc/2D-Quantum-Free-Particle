"""Deterministic restartable Adam runner for constrained inverse design."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any
import json
import platform
import subprocess

import jax
import jax.numpy as jnp
import numpy as np

from . import __version__
from .design_config import DesignModeConfig, InverseDesignConfig
from .jax_constraints import parameterize_phase_jax, phase_total_variation_jax
from .jax_objectives import (
    intensity_similarity_jax,
    mode_overlap_jax,
    mode_power_efficiency_jax,
)
from .jax_propagation import multiplane_scan_jax
from .optics import gaussian_mode, hermite_gaussian_mode, make_optical_grid
from .optics_propagation import fresnel_transfer_function


@dataclass
class _OptimizationState:
    parameters: Any
    first_moment: Any
    second_moment: Any
    best_parameters: Any
    step: int
    best_loss: float
    best_iteration: int
    patience_loss: float
    stall_count: int
    history: list[dict[str, Any]]


@dataclass(frozen=True)
class OptimizationResult:
    out_dir: Path
    status: str
    metrics: dict[str, Any]


def run_inverse_design(
    config: InverseDesignConfig,
    out_dir: str | Path,
    *,
    resume: str | Path | None = None,
    max_steps: int | None = None,
) -> OptimizationResult:
    """Run or resume deterministic constrained Adam optimization."""

    jax.config.update("jax_enable_x64", True)
    _validate_config(config)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fingerprint = _config_fingerprint(config)
    grid = make_optical_grid(
        config.grid.shape,
        config.grid.extent_m,
        config.grid.wavelength_m,
    )
    source = _make_mode(config.input_mode, grid)
    target = _make_mode(config.target_mode, grid)
    transfers = np.stack(
        [
            fresnel_transfer_function(grid, distance)
            for distance in config.design.plane_distances_m
        ]
    )
    source_jax = jnp.asarray(source)
    target_jax = jnp.asarray(target)
    transfers_jax = jnp.asarray(transfers)
    area = grid.sample_area

    def evaluate(parameters, *, straight_through: bool):
        phases = parameterize_phase_jax(
            parameters,
            bounds=config.design.phase_bounds_rad,
            smoothing_passes=config.design.smoothing_passes,
            quantization_levels=config.design.quantization_levels,
            straight_through=straight_through,
        )
        final_field, _ = multiplane_scan_jax(source_jax, phases, transfers_jax)
        overlap = mode_overlap_jax(final_field, target_jax, area)
        intensity = intensity_similarity_jax(final_field, target_jax, area)
        efficiency = mode_power_efficiency_jax(final_field, target_jax, 1.0, area)
        total_variation = jnp.mean(
            jax.vmap(phase_total_variation_jax)(phases)
        )
        loss = (
            config.objective.mode_overlap_weight * (1.0 - overlap)
            + config.objective.intensity_weight * (1.0 - intensity)
            + config.design.total_variation_weight * total_variation
        )
        auxiliary = jnp.asarray(
            [overlap, intensity, efficiency, total_variation],
            dtype=final_field.real.dtype,
        )
        return loss, auxiliary

    def training_objective(parameters):
        return evaluate(parameters, straight_through=True)

    def hard_objective(parameters):
        return evaluate(parameters, straight_through=False)

    @jax.jit
    def adam_step(parameters, first_moment, second_moment, step):
        (_, _), gradient = jax.value_and_grad(
            training_objective,
            has_aux=True,
        )(parameters)
        first_moment = (
            config.optimizer.beta1 * first_moment
            + (1.0 - config.optimizer.beta1) * gradient
        )
        second_moment = (
            config.optimizer.beta2 * second_moment
            + (1.0 - config.optimizer.beta2) * gradient**2
        )
        first_unbiased = first_moment / (1.0 - config.optimizer.beta1**step)
        second_unbiased = second_moment / (1.0 - config.optimizer.beta2**step)
        parameters = parameters - config.optimizer.learning_rate * first_unbiased / (
            jnp.sqrt(second_unbiased) + config.optimizer.epsilon
        )
        training_loss, _ = training_objective(parameters)
        return (
            parameters,
            first_moment,
            second_moment,
            training_loss,
            jnp.linalg.norm(gradient),
        )

    hard_evaluate = jax.jit(hard_objective)

    def fields(parameters):
        phases = parameterize_phase_jax(
            parameters,
            bounds=config.design.phase_bounds_rad,
            smoothing_passes=config.design.smoothing_passes,
            quantization_levels=config.design.quantization_levels,
            straight_through=False,
        )
        final_field, history = multiplane_scan_jax(
            source_jax,
            phases,
            transfers_jax,
        )
        return phases, final_field, history

    compiled_fields = jax.jit(fields)
    if resume is None:
        rng = np.random.default_rng(config.seed)
        parameters = jnp.asarray(
            config.design.initial_scale
            * rng.standard_normal((len(transfers), *grid.shape)),
            dtype=jnp.float64,
        )
        initial_loss, initial_aux = hard_evaluate(parameters)
        initial_record = _history_record(
            0,
            float(initial_loss),
            initial_aux,
            float(initial_loss),
            0.0,
        )
        state = _OptimizationState(
            parameters=parameters,
            first_moment=jnp.zeros_like(parameters),
            second_moment=jnp.zeros_like(parameters),
            best_parameters=parameters,
            step=0,
            best_loss=float(initial_loss),
            best_iteration=0,
            patience_loss=float(initial_loss),
            stall_count=0,
            history=[initial_record],
        )
        _save_best_design(out_dir, state, compiled_fields)
    else:
        state = _load_checkpoint(Path(resume), fingerprint)

    target_step = config.optimizer.iterations
    if max_steps is not None:
        if max_steps < 0:
            raise ValueError("max_steps must be nonnegative")
        target_step = min(target_step, state.step + max_steps)

    early_stopped = False
    while state.step < target_step:
        iteration = state.step + 1
        (
            state.parameters,
            state.first_moment,
            state.second_moment,
            training_loss,
            gradient_norm,
        ) = adam_step(
            state.parameters,
            state.first_moment,
            state.second_moment,
            jnp.asarray(iteration, dtype=jnp.float64),
        )
        hard_loss, hard_aux = hard_evaluate(state.parameters)
        loss_value = float(hard_loss)
        state.step = iteration
        state.history.append(
            _history_record(
                iteration,
                loss_value,
                hard_aux,
                float(training_loss),
                float(gradient_norm),
            )
        )
        if loss_value < state.best_loss:
            state.best_loss = loss_value
            state.best_iteration = iteration
            state.best_parameters = state.parameters
            _save_best_design(out_dir, state, compiled_fields)
        if loss_value < state.patience_loss - config.optimizer.min_delta:
            state.patience_loss = loss_value
            state.stall_count = 0
        else:
            state.stall_count += 1
        if iteration % config.optimizer.checkpoint_interval == 0:
            _save_checkpoint(out_dir / "checkpoint.npz", state, fingerprint)
        if (
            config.optimizer.early_stopping_patience > 0
            and state.stall_count >= config.optimizer.early_stopping_patience
        ):
            early_stopped = True
            break

    _save_checkpoint(out_dir / "checkpoint.npz", state, fingerprint)
    if early_stopped:
        status = "early_stopped"
    elif state.step >= config.optimizer.iterations:
        status = "complete"
    else:
        status = "incomplete"
    metrics = _write_artifacts(
        config,
        out_dir,
        status,
        state,
        source,
        target,
        grid,
        compiled_fields,
        hard_evaluate,
        fingerprint,
    )
    return OptimizationResult(out_dir=out_dir, status=status, metrics=metrics)


def _history_record(
    iteration: int,
    loss: float,
    auxiliary: Any,
    training_loss: float,
    gradient_norm: float,
) -> dict[str, Any]:
    overlap, intensity, efficiency, total_variation = np.asarray(auxiliary)
    return {
        "iteration": iteration,
        "loss": loss,
        "training_loss": training_loss,
        "mode_overlap": float(overlap),
        "intensity_similarity": float(intensity),
        "mode_power_efficiency": float(efficiency),
        "total_variation": float(total_variation),
        "gradient_norm": gradient_norm,
    }


def _make_mode(mode: DesignModeConfig, grid) -> np.ndarray:
    kind = mode.kind.lower()
    if kind == "gaussian":
        return gaussian_mode(grid, mode.waist_m, center=mode.center_m)
    if kind in {"hermite_gaussian", "hg"}:
        return hermite_gaussian_mode(
            grid,
            mode.order,
            mode.waist_m,
            center=mode.center_m,
        )
    raise ValueError(f"unsupported design mode kind: {mode.kind}")


def _save_checkpoint(
    path: Path,
    state: _OptimizationState,
    fingerprint: str,
) -> None:
    np.savez_compressed(
        path,
        config_fingerprint=fingerprint,
        parameters=np.asarray(state.parameters),
        first_moment=np.asarray(state.first_moment),
        second_moment=np.asarray(state.second_moment),
        best_parameters=np.asarray(state.best_parameters),
        step=state.step,
        best_loss=state.best_loss,
        best_iteration=state.best_iteration,
        patience_loss=state.patience_loss,
        stall_count=state.stall_count,
        history=json.dumps(state.history, sort_keys=True),
    )


def _load_checkpoint(path: Path, fingerprint: str) -> _OptimizationState:
    with np.load(path) as data:
        saved_fingerprint = str(data["config_fingerprint"])
        if saved_fingerprint != fingerprint:
            raise ValueError("checkpoint config fingerprint does not match config")
        return _OptimizationState(
            parameters=jnp.asarray(data["parameters"]),
            first_moment=jnp.asarray(data["first_moment"]),
            second_moment=jnp.asarray(data["second_moment"]),
            best_parameters=jnp.asarray(data["best_parameters"]),
            step=int(data["step"]),
            best_loss=float(data["best_loss"]),
            best_iteration=int(data["best_iteration"]),
            patience_loss=float(data["patience_loss"]),
            stall_count=int(data["stall_count"]),
            history=json.loads(str(data["history"])),
        )


def _save_best_design(out_dir: Path, state: _OptimizationState, compiled_fields) -> None:
    phases, output, plane_fields = compiled_fields(state.best_parameters)
    np.savez_compressed(
        out_dir / "best_design.npz",
        parameters=np.asarray(state.best_parameters),
        phase_masks=np.asarray(phases),
        output_field=np.asarray(output),
        plane_fields=np.asarray(plane_fields),
        iteration=state.best_iteration,
        loss=state.best_loss,
    )


def _write_artifacts(
    config: InverseDesignConfig,
    out_dir: Path,
    status: str,
    state: _OptimizationState,
    source: np.ndarray,
    target: np.ndarray,
    grid,
    compiled_fields,
    hard_evaluate,
    fingerprint: str,
) -> dict[str, Any]:
    best_phases, best_output, best_plane_fields = compiled_fields(state.best_parameters)
    final_phases, final_output, final_plane_fields = compiled_fields(state.parameters)
    best_loss, best_aux = hard_evaluate(state.best_parameters)
    best_record = _history_record(
        state.best_iteration,
        float(best_loss),
        best_aux,
        float(best_loss),
        0.0,
    )
    initial = state.history[0]
    gradient_evidence = _gradient_evidence()
    metrics = {
        "name": config.name,
        "status": status,
        "seed": config.seed,
        "backend": "jax",
        "device": jax.default_backend(),
        "jax_x64": bool(jax.config.x64_enabled),
        "config_fingerprint": fingerprint,
        "iterations_completed": state.step,
        "best_iteration": state.best_iteration,
        "initial": initial,
        "best": best_record,
        "final": state.history[-1],
        "mode_overlap_improvement": best_record["mode_overlap"]
        - initial["mode_overlap"],
        "loss_reduction": initial["loss"] - best_record["loss"],
        "constraints": config.to_dict()["design"],
        "objective": config.to_dict()["objective"],
        "gradient_evidence": gradient_evidence,
    }
    environment = {
        "created_at": datetime.now(UTC).isoformat(),
        "package_version": __version__,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "jax_version": jax.__version__,
        "dependencies": {
            name: _package_version(name)
            for name in ("jaxlib", "matplotlib", "numpy", "pytest", "scipy")
        },
        "device": str(jax.devices()[0]),
        "git_commit": _git_value("rev-parse", "HEAD"),
        "git_dirty": bool(_git_value("status", "--porcelain")),
        "seed": config.seed,
    }
    (out_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "history.json").write_text(
        json.dumps(state.history, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "resolved_config.json").write_text(
        json.dumps(config.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    np.savez_compressed(
        out_dir / "fields.npz",
        x=grid.x,
        y=grid.y,
        input_field=source,
        target_field=target,
        best_output_field=np.asarray(best_output),
        final_output_field=np.asarray(final_output),
        best_phase_masks=np.asarray(best_phases),
        final_phase_masks=np.asarray(final_phases),
        best_plane_fields=np.asarray(best_plane_fields),
        final_plane_fields=np.asarray(final_plane_fields),
    )
    (out_dir / "report.md").write_text(
        _report_markdown(config, metrics, environment),
        encoding="utf-8",
    )
    return metrics


def _gradient_evidence() -> dict[str, Any]:
    path = Path(__file__).resolve().parents[2] / "benchmarks" / "reference" / "jax_m2.json"
    if not path.exists():
        return {"path": str(path), "available": False}
    evidence = json.loads(path.read_text(encoding="utf-8"))
    return {
        "path": "benchmarks/reference/jax_m2.json",
        "available": True,
        "seed": evidence["seed"],
        "declared_limits": evidence["declared_limits"],
        "mode_relative_error": evidence["metrics"]["mode_gradient"][
            "relative_error"
        ],
        "intensity_relative_error": evidence["metrics"]["intensity_gradient"][
            "relative_error"
        ],
    }


def _report_markdown(
    config: InverseDesignConfig,
    metrics: dict[str, Any],
    environment: dict[str, Any],
) -> str:
    return f"""# {config.name} inverse-design report

Status: **{metrics['status']}** after {metrics['iterations_completed']} updates.

## Objective

- Initial mode overlap: {metrics['initial']['mode_overlap']:.8f}
- Best mode overlap: {metrics['best']['mode_overlap']:.8f}
- Mode-overlap improvement: {metrics['mode_overlap_improvement']:.8f}
- Best intensity similarity: {metrics['best']['intensity_similarity']:.8f}
- Best input-referenced mode efficiency: {metrics['best']['mode_power_efficiency']:.8f}
- Total loss reduction: {metrics['loss_reduction']:.8f}

## Constraints and reproducibility

- Seed: {config.seed}
- Phase bounds: {config.design.phase_bounds_rad}
- Smoothing passes: {config.design.smoothing_passes}
- TV weight: {config.design.total_variation_weight}
- Quantization levels: {config.design.quantization_levels}
- Backend/device: JAX x64 / {metrics['device']}
- Git commit: {environment['git_commit']} (dirty={environment['git_dirty']})
- Config fingerprint: `{metrics['config_fingerprint']}`

M2 centered directional gradient evidence is recorded in
`benchmarks/reference/jax_m2.json`; the report metrics include its declared
tolerances and measured mode/intensity relative errors.

## Scientific scope

This optimization uses a monochromatic coherent scalar Fresnel model. Objective
improvement validates the deterministic discrete optimization workflow only. It
does not establish robustness, band-limited angular-spectrum transfer, vector-
Maxwell accuracy, or manufacturing readiness.
"""


def _config_fingerprint(config: InverseDesignConfig) -> str:
    payload = json.dumps(config.to_dict(), sort_keys=True, separators=(",", ":"))
    return sha256(payload.encode("utf-8")).hexdigest()


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


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _validate_config(config: InverseDesignConfig) -> None:
    if not config.design.plane_distances_m or any(
        distance <= 0.0 for distance in config.design.plane_distances_m
    ):
        raise ValueError("plane distances must be nonempty and positive")
    if not np.isfinite(config.design.initial_scale) or config.design.initial_scale < 0.0:
        raise ValueError("initial_scale must be finite and nonnegative")
    lower, upper = config.design.phase_bounds_rad
    if lower >= upper:
        raise ValueError("phase bounds must be strictly increasing")
    if config.design.smoothing_passes < 0:
        raise ValueError("smoothing_passes must be nonnegative")
    if config.design.total_variation_weight < 0.0:
        raise ValueError("total_variation_weight must be nonnegative")
    if (
        config.design.quantization_levels is not None
        and config.design.quantization_levels < 2
    ):
        raise ValueError("quantization_levels must be at least two")
    if config.optimizer.iterations < 1:
        raise ValueError("optimizer iterations must be positive")
    if config.optimizer.learning_rate <= 0.0:
        raise ValueError("learning rate must be positive")
    if not 0.0 < config.optimizer.beta1 < 1.0:
        raise ValueError("Adam beta1 must lie between zero and one")
    if not 0.0 < config.optimizer.beta2 < 1.0:
        raise ValueError("Adam beta2 must lie between zero and one")
    if config.optimizer.epsilon <= 0.0:
        raise ValueError("Adam epsilon must be positive")
    if config.optimizer.checkpoint_interval < 1:
        raise ValueError("checkpoint_interval must be positive")
    if config.optimizer.early_stopping_patience < 0:
        raise ValueError("early_stopping_patience must be nonnegative")
    if config.optimizer.min_delta < 0.0:
        raise ValueError("min_delta must be nonnegative")
    if (
        config.objective.mode_overlap_weight < 0.0
        or config.objective.intensity_weight < 0.0
    ):
        raise ValueError("objective weights must be nonnegative")
    if (
        config.objective.mode_overlap_weight == 0.0
        and config.objective.intensity_weight == 0.0
        and config.design.total_variation_weight == 0.0
    ):
        raise ValueError("at least one objective or regularization weight must be positive")
