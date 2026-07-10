from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any
import json
import subprocess

import numpy as np

from . import __version__
from .experiments import ExperimentResult


def save_run(result: ExperimentResult, out_dir: str | Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_path = out_dir / "run.npz"
    metadata = _metadata(result)
    payload: dict[str, Any] = {
        "x": result.grid.x,
        "y": result.grid.y,
        "potential": result.potential,
        "metadata": json.dumps(metadata, indent=2, sort_keys=True),
    }
    payload.update(result.arrays)
    for name, mask in result.masks.items():
        payload[f"mask_{name}"] = mask
    np.savez_compressed(run_path, **payload)
    (out_dir / "metrics.json").write_text(
        json.dumps(result.metrics, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return run_path


def load_run(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with np.load(path) as data:
        loaded: dict[str, Any] = {key: data[key] for key in data.files}
    loaded["metadata"] = json.loads(str(loaded["metadata"]))
    return loaded


def write_metrics(path: str | Path, metrics: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8")


def _metadata(result: ExperimentResult) -> dict[str, Any]:
    return {
        "backend": str(result.metrics.get("backend", "numpy")),
        "config": result.config.to_dict(),
        "created_at": datetime.now(UTC).isoformat(),
        "git_commit": _git_commit(),
        "package_version": __version__,
        "units": {"hbar": result.config.solver.hbar, "mass": result.config.solver.mass},
    }


def _git_commit() -> str | None:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()
