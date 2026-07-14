"""Regenerate deterministic M2 objective and gradient evidence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np

from quantum_dynamics_lab.jax_objectives import (
    intensity_similarity_jax,
    mode_overlap_jax,
    mode_power_efficiency_jax,
)
from quantum_dynamics_lab.jax_propagation import multiplane_scan_jax
from quantum_dynamics_lab.objectives import (
    intensity_similarity,
    mode_overlap,
    mode_power_efficiency,
)
from quantum_dynamics_lab.optics import (
    gaussian_mode,
    hermite_gaussian_mode,
    make_optical_grid,
)
from quantum_dynamics_lab.optics_propagation import fresnel_transfer_function


REFERENCE_PATH = Path(__file__).with_name("jax_m2.json")
SEED = 20260714
FINITE_DIFFERENCE_STEP = 1.0e-5


def _directional_check(
    loss: Callable[[jax.Array], jax.Array],
    parameters: jax.Array,
    direction: jax.Array,
) -> dict[str, float]:
    value, gradient = jax.value_and_grad(loss)(parameters)
    autodiff = jnp.vdot(gradient, direction).real
    plus = loss(parameters + FINITE_DIFFERENCE_STEP * direction)
    minus = loss(parameters - FINITE_DIFFERENCE_STEP * direction)
    finite_difference = (plus - minus) / (2.0 * FINITE_DIFFERENCE_STEP)
    absolute_error = jnp.abs(autodiff - finite_difference)
    scale = jnp.maximum(
        jnp.maximum(jnp.abs(autodiff), jnp.abs(finite_difference)),
        1.0e-12,
    )
    return {
        "loss": float(value),
        "autodiff_directional_derivative": float(autodiff),
        "finite_difference_directional_derivative": float(finite_difference),
        "absolute_error": float(absolute_error),
        "relative_error": float(absolute_error / scale),
    }


def generate_metrics() -> dict[str, Any]:
    grid = make_optical_grid(32, 4.0e-3, 1064.0e-9)
    source = gaussian_mode(grid, 0.55e-3)
    target = hermite_gaussian_mode(grid, (1, 0), 0.55e-3)
    distances = (0.04, 0.05, 0.03)
    transfers = np.stack(
        [fresnel_transfer_function(grid, distance) for distance in distances]
    )
    rng = np.random.default_rng(SEED)
    parameters = 0.25 * rng.standard_normal((len(distances), *grid.shape))
    direction = rng.standard_normal(parameters.shape)
    direction /= np.linalg.norm(direction)

    source_jax = jnp.asarray(source)
    target_jax = jnp.asarray(target)
    transfers_jax = jnp.asarray(transfers)
    parameters_jax = jnp.asarray(parameters)
    direction_jax = jnp.asarray(direction)

    def final_field(phase_masks):
        return multiplane_scan_jax(
            source_jax,
            phase_masks,
            transfers_jax,
        )[0]

    def mode_loss(phase_masks):
        return 1.0 - mode_overlap_jax(
            final_field(phase_masks),
            target_jax,
            grid.sample_area,
        )

    def intensity_loss(phase_masks):
        return 1.0 - intensity_similarity_jax(
            final_field(phase_masks),
            target_jax,
            grid.sample_area,
        )

    final = final_field(parameters_jax)
    final_numpy = np.asarray(final)
    numpy_values = {
        "mode_overlap": mode_overlap(final_numpy, target, grid.sample_area),
        "intensity_similarity": intensity_similarity(
            final_numpy,
            target,
            grid.sample_area,
        ),
        "mode_power_efficiency": mode_power_efficiency(
            final_numpy,
            target,
            1.0,
            grid.sample_area,
        ),
    }
    jax_values = {
        "mode_overlap": float(mode_overlap_jax(final, target_jax, grid.sample_area)),
        "intensity_similarity": float(
            intensity_similarity_jax(final, target_jax, grid.sample_area)
        ),
        "mode_power_efficiency": float(
            mode_power_efficiency_jax(final, target_jax, 1.0, grid.sample_area)
        ),
    }

    @jax.jit
    def compiled_update(phase_masks):
        value, gradient = jax.value_and_grad(mode_loss)(phase_masks)
        return phase_masks - 0.05 * gradient, value, jnp.linalg.norm(gradient)

    updated, update_loss, gradient_norm = compiled_update(parameters_jax)
    updated.block_until_ready()
    return {
        "schema_version": 1,
        "seed": SEED,
        "configuration": {
            "grid": 32,
            "extent_m": 4.0e-3,
            "wavelength_m": 1064.0e-9,
            "waist_m": 0.55e-3,
            "planes": len(distances),
            "distances_m": list(distances),
            "finite_difference_step": FINITE_DIFFERENCE_STEP,
            "jax_platform": jax.default_backend(),
            "jax_x64": bool(jax.config.x64_enabled),
        },
        "metrics": {
            "objective_value_max_abs_error": max(
                abs(numpy_values[name] - jax_values[name]) for name in numpy_values
            ),
            "mode_gradient": _directional_check(
                mode_loss,
                parameters_jax,
                direction_jax,
            ),
            "intensity_gradient": _directional_check(
                intensity_loss,
                parameters_jax,
                direction_jax,
            ),
            "compiled_update_loss": float(update_loss),
            "compiled_update_gradient_norm": float(gradient_norm),
            "compiled_update_all_finite": bool(np.all(np.isfinite(np.asarray(updated)))),
        },
        "declared_limits": {
            "objective_value_max_abs_error": 1.0e-12,
            "gradient_absolute_error": 1.0e-9,
            "gradient_relative_error": 1.0e-6,
        },
    }


def main() -> None:
    REFERENCE_PATH.write_text(
        json.dumps(generate_metrics(), indent=2) + "\n",
        encoding="utf-8",
    )
    print(REFERENCE_PATH)


if __name__ == "__main__":
    main()
