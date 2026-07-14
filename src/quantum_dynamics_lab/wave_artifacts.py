"""Artifact persistence for scalar-optics forward runs."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
import json
import platform
import subprocess

import numpy as np

from . import __version__
from .wave_experiments import WaveRun


def save_wave_run(run: WaveRun, out_dir: str | Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_path = out_dir / "wave_run.npz"
    metadata = {
        "created_at": datetime.now(UTC).isoformat(),
        "package_version": __version__,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "git_commit": _git_value("rev-parse", "HEAD"),
        "git_dirty": bool(_git_value("status", "--porcelain")),
        "config": run.config.to_dict(),
    }
    np.savez_compressed(
        run_path,
        x=run.grid.x,
        y=run.grid.y,
        input_field=run.input_field,
        final_field=run.final_field,
        metadata=json.dumps(metadata, indent=2, sort_keys=True),
    )
    (out_dir / "metrics.json").write_text(
        json.dumps(run.metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "resolved_config.json").write_text(
        json.dumps(run.config.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return run_path


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
