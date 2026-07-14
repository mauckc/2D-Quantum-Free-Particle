from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest


jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)

from quantum_dynamics_lab.design_config import load_inverse_design_config
from quantum_dynamics_lab.optimization import run_inverse_design
from quantum_dynamics_lab.wave_cli import main as wave_main


REPOSITORY_ROOT = Path(__file__).parents[1]
SMOKE_CONFIG = REPOSITORY_ROOT / "examples" / "inverse_design_smoke.toml"


def test_smoke_optimization_improves_and_writes_complete_artifacts(
    tmp_path: Path,
) -> None:
    config = load_inverse_design_config(SMOKE_CONFIG)
    result = run_inverse_design(config, tmp_path)

    assert result.status in {"complete", "early_stopped"}
    assert result.metrics["mode_overlap_improvement"] > 0.2
    assert result.metrics["loss_reduction"] > 0.2
    assert result.metrics["best"]["mode_overlap"] > result.metrics["initial"][
        "mode_overlap"
    ]
    for name in (
        "best_design.npz",
        "checkpoint.npz",
        "environment.json",
        "fields.npz",
        "history.json",
        "metrics.json",
        "report.md",
        "resolved_config.json",
    ):
        assert (tmp_path / name).exists()
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "M2 centered directional gradient evidence" in report
    assert "monochromatic coherent scalar Fresnel" in report
    with np.load(tmp_path / "best_design.npz") as design:
        assert design["phase_masks"].shape == (
            len(config.design.plane_distances_m),
            *config.grid.shape,
        )


def test_interrupted_resume_matches_uninterrupted_history(tmp_path: Path) -> None:
    base = load_inverse_design_config(SMOKE_CONFIG)
    config = replace(
        base,
        optimizer=replace(
            base.optimizer,
            iterations=8,
            checkpoint_interval=4,
            early_stopping_patience=0,
        ),
    )
    full = run_inverse_design(config, tmp_path / "full")
    partial = run_inverse_design(config, tmp_path / "resumed", max_steps=4)
    assert partial.status == "incomplete"
    resumed = run_inverse_design(
        config,
        tmp_path / "resumed",
        resume=tmp_path / "resumed" / "checkpoint.npz",
    )

    assert full.status == resumed.status == "complete"
    full_history = json.loads(
        (tmp_path / "full" / "history.json").read_text(encoding="utf-8")
    )
    resumed_history = json.loads(
        (tmp_path / "resumed" / "history.json").read_text(encoding="utf-8")
    )
    assert resumed_history == full_history
    assert resumed.metrics["best"] == full.metrics["best"]
    with np.load(tmp_path / "full" / "checkpoint.npz") as full_checkpoint:
        with np.load(
            tmp_path / "resumed" / "checkpoint.npz"
        ) as resumed_checkpoint:
            np.testing.assert_array_equal(
                resumed_checkpoint["parameters"],
                full_checkpoint["parameters"],
            )


def test_early_stopping_fingerprint_guard_and_cli(tmp_path: Path) -> None:
    base = load_inverse_design_config(SMOKE_CONFIG)
    early_config = replace(
        base,
        optimizer=replace(
            base.optimizer,
            iterations=20,
            early_stopping_patience=2,
            min_delta=10.0,
            checkpoint_interval=1,
        ),
    )
    early = run_inverse_design(early_config, tmp_path / "early")
    assert early.status == "early_stopped"
    assert early.metrics["iterations_completed"] == 2
    changed = replace(early_config, seed=early_config.seed + 1)
    with pytest.raises(ValueError, match="fingerprint"):
        run_inverse_design(
            changed,
            tmp_path / "early",
            resume=tmp_path / "early" / "checkpoint.npz",
        )

    cli_config = replace(
        base,
        optimizer=replace(
            base.optimizer,
            iterations=2,
            early_stopping_patience=0,
            checkpoint_interval=1,
        ),
    )
    config_path = tmp_path / "cli.toml"
    config_path.write_text(
        SMOKE_CONFIG.read_text(encoding="utf-8")
        .replace("iterations = 30", "iterations = 2")
        .replace("early_stopping_patience = 15", "early_stopping_patience = 0")
        .replace("checkpoint_interval = 5", "checkpoint_interval = 1"),
        encoding="utf-8",
    )
    assert cli_config.optimizer.iterations == 2
    assert wave_main(
        ["optimize", str(config_path), "--out", str(tmp_path / "cli")]
    ) == 0


def test_saved_quantized_best_design_uses_hard_levels(tmp_path: Path) -> None:
    base = load_inverse_design_config(SMOKE_CONFIG)
    config = replace(
        base,
        design=replace(base.design, quantization_levels=4),
        optimizer=replace(
            base.optimizer,
            iterations=2,
            checkpoint_interval=1,
            early_stopping_patience=0,
        ),
    )

    result = run_inverse_design(config, tmp_path)

    assert result.status == "complete"
    with np.load(tmp_path / "best_design.npz") as design:
        phase_masks = design["phase_masks"]
    expected = -np.pi + np.arange(4) * (2.0 * np.pi / 4)
    distances = np.min(
        np.abs(np.unique(phase_masks)[:, None] - expected[None, :]),
        axis=1,
    )
    np.testing.assert_allclose(distances, 0.0, atol=1e-15)
