from __future__ import annotations

from pathlib import Path
from typing import Any
import shutil

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np


def render_run(run_path: str | Path, out_dir: str | Path | None = None) -> dict[str, Path]:
    from .artifacts import load_run

    run_path = Path(run_path)
    out_dir = Path(out_dir) if out_dir is not None else run_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    data = load_run(run_path)
    outputs: dict[str, Path] = {}
    outputs["potential"] = _plot_potential(data, out_dir / "potential.png")
    outputs["summary"] = _plot_summary(data, out_dir / "summary.png")
    if data["metadata"]["config"]["output"].get("render_video", True):
        video = _render_video(data, out_dir / "probability.mp4")
        if video is not None:
            outputs["video"] = video
    return outputs


def render_comparison(run_dirs: list[Path], out_dir: str | Path) -> Path:
    from .artifacts import load_run

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for run_dir in run_dirs:
        run = load_run(run_dir / "run.npz")
        metadata = run["metadata"]
        config = metadata["config"]
        interval = config["measurement"].get("interval")
        detected = _last(run.get("detected_probability"))
        transmitted = _last(run.get("transmission_probability"))
        rows.append((config["name"], interval or 0.0, detected, transmitted))
    rows.sort(key=lambda row: row[1])

    fig, ax = plt.subplots(figsize=(7, 4))
    labels = [row[0] for row in rows]
    detected = [row[2] for row in rows]
    transmitted = [row[3] for row in rows]
    x = np.arange(len(rows))
    ax.bar(x - 0.18, detected, width=0.36, label="Detected")
    ax.bar(x + 0.18, transmitted, width=0.36, label="Conditional transmitted")
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("Probability")
    ax.set_title("Zeno measurement sweep")
    ax.legend()
    fig.tight_layout()
    path = out_dir / "comparison.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def _plot_potential(data: dict[str, Any], path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(
        data["potential"],
        origin="lower",
        extent=_extent(data),
        cmap="magma",
        aspect="equal",
    )
    ax.set_title("Potential")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(image, ax=ax, label="V(x, y)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def _plot_summary(data: dict[str, Any], path: Path) -> Path:
    psi_frames = data["psi_frames"]
    final_probability = np.abs(psi_frames[-1]) ** 2
    final_phase = np.angle(psi_frames[-1])
    metric_times = data.get("metric_times", data["times"])
    norm = data.get("norm")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    prob_image = axes[0, 0].imshow(
        final_probability,
        origin="lower",
        extent=_extent(data),
        cmap="viridis",
        aspect="equal",
    )
    axes[0, 0].set_title("Final probability density")
    fig.colorbar(prob_image, ax=axes[0, 0])

    phase_image = axes[0, 1].imshow(
        final_phase,
        origin="lower",
        extent=_extent(data),
        cmap="twilight",
        aspect="equal",
    )
    axes[0, 1].set_title("Final phase")
    fig.colorbar(phase_image, ax=axes[0, 1])

    if norm is not None:
        axes[1, 0].plot(metric_times, norm, label="coherent")
    if "which_path_norm" in data:
        axes[1, 0].plot(metric_times, data["which_path_norm"], label="which-path")
    axes[1, 0].set_title("Norm")
    axes[1, 0].set_xlabel("time")
    axes[1, 0].set_ylabel("probability mass")
    axes[1, 0].legend(loc="best")

    if "coherent_screen_profile" in data:
        y = data["y"]
        axes[1, 1].plot(y, data["coherent_screen_profile"], label="coherent")
        axes[1, 1].plot(y, data["which_path_screen_profile"], label="which-path")
        axes[1, 1].set_title("Screen profile")
        axes[1, 1].set_xlabel("y")
        axes[1, 1].legend(loc="best")
    else:
        if "transmission_probability" in data:
            axes[1, 1].plot(
                metric_times,
                data["transmission_probability"],
                label="conditional transmitted",
            )
        if "detected_probability" in data:
            axes[1, 1].plot(metric_times, data["detected_probability"], label="detected")
        if "survival_weight" in data:
            axes[1, 1].plot(metric_times, data["survival_weight"], label="survival")
        axes[1, 1].set_title("Measurement probabilities")
        axes[1, 1].set_xlabel("time")
        axes[1, 1].legend(loc="best")

    for ax in axes.flat[:2]:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def _render_video(data: dict[str, Any], path: Path) -> Path | None:
    if shutil.which("ffmpeg") is None:
        return None
    frames = np.abs(data["psi_frames"]) ** 2
    vmax = float(np.max(frames))
    if vmax <= 0:
        vmax = 1.0
    fps = int(data["metadata"]["config"]["output"].get("fps", 20))
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(
        frames[0],
        origin="lower",
        extent=_extent(data),
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        aspect="equal",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = ax.set_title("t = 0.000")
    fig.colorbar(image, ax=ax, label="|psi|^2")
    writer = FFMpegWriter(fps=fps)
    with writer.saving(fig, str(path), dpi=120):
        for time, frame in zip(data["times"], frames):
            image.set_data(frame)
            title.set_text(f"t = {float(time):.3f}")
            writer.grab_frame()
    plt.close(fig)
    return path


def _extent(data: dict[str, Any]) -> tuple[float, float, float, float]:
    x = data["x"]
    y = data["y"]
    dx = float(x[1] - x[0]) if len(x) > 1 else 1.0
    dy = float(y[1] - y[0]) if len(y) > 1 else 1.0
    return (
        float(x[0] - dx / 2),
        float(x[-1] + dx / 2),
        float(y[0] - dy / 2),
        float(y[-1] + dy / 2),
    )


def _last(value: Any) -> float:
    if value is None:
        return 0.0
    return float(np.asarray(value)[-1])
