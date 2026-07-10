from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
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
    outputs["report"] = _write_run_report(data, outputs, out_dir / "report.md")
    if data["metadata"]["config"]["output"].get("render_video", True):
        video = _render_video(data, out_dir / "probability.mp4")
        if video is not None:
            outputs["video"] = video
    return outputs


def render_comparison(run_dirs: list[Path], out_dir: str | Path) -> Path:
    from .artifacts import load_run

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | str | None]] = []
    for run_dir in run_dirs:
        run = load_run(run_dir / "run.npz")
        metadata = run["metadata"]
        config = metadata["config"]
        rows.append(_comparison_row(config, run))
    rows.sort(key=lambda row: (float(row["barrier_height"]), float(row["interval_sort"])))
    _write_comparison_csv(rows, out_dir / "comparison.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=False)
    heights = sorted({float(row["barrier_height"]) for row in rows})
    for height in heights:
        group = [row for row in rows if float(row["barrier_height"]) == height]
        group.sort(key=lambda row: float(row["interval_sort"]))
        x = [float(row["interval_sort"]) for row in group]
        label = f"V0={height:g}"
        axes[0].plot(
            x,
            [float(row["unconditional_transmission_probability"]) for row in group],
            marker="o",
            label=label,
        )
        axes[1].plot(
            x,
            [float(row["survival_weight"]) for row in group],
            marker="o",
            label=label,
        )
    axes[0].set_title("No-click branch transmission")
    axes[0].set_xlabel("measurement interval (0 = unmeasured)")
    axes[0].set_ylabel("probability")
    axes[1].set_title("No-click survival weight")
    axes[1].set_xlabel("measurement interval (0 = unmeasured)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend()
    fig.tight_layout()
    path = out_dir / "comparison.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    heatmap_path = _plot_transmission_heatmap(rows, out_dir / "zeno_transmission_heatmap.png")
    _write_comparison_report(rows, path, heatmap_path, out_dir / "report.md")
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
        if "unconditional_transmission_probability" in data:
            axes[1, 1].plot(
                metric_times,
                data["unconditional_transmission_probability"],
                label="unconditional transmitted",
            )
        if "detected_probability" in data:
            axes[1, 1].plot(metric_times, data["detected_probability"], label="detected")
        if "absorbed_probability" in data:
            axes[1, 1].plot(metric_times, data["absorbed_probability"], label="absorbed")
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


def _comparison_row(config: dict[str, Any], run: dict[str, Any]) -> dict[str, float | str | None]:
    interval = config["measurement"].get("interval")
    interval_sort = 0.0 if interval is None else float(interval)
    return {
        "name": config["name"],
        "barrier_height": float(config["potential"]["height"]),
        "interval": None if interval is None else float(interval),
        "interval_sort": interval_sort,
        "detector_x": _detector_x_from_config(config),
        "detected_probability": _last(run.get("detected_probability")),
        "absorbed_probability": _last(run.get("absorbed_probability")),
        "survival_weight": _last(run.get("survival_weight")),
        "conditional_transmission_probability": _last(run.get("transmission_probability")),
        "unconditional_transmission_probability": _last(
            run.get("unconditional_transmission_probability")
        ),
    }


def _write_comparison_csv(rows: list[dict[str, float | str | None]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _detector_x_from_config(config: dict[str, Any]) -> float:
    measurement = config["measurement"]
    if measurement.get("detector_x") is not None:
        return float(measurement["detector_x"])
    potential = config["potential"]
    return (
        float(potential["barrier_x"])
        + float(potential["barrier_width"]) / 2.0
        + float(measurement.get("detector_buffer", 0.0))
    )


def _plot_transmission_heatmap(rows: list[dict[str, float | str | None]], path: Path) -> Path:
    heights = sorted({float(row["barrier_height"]) for row in rows})
    intervals = sorted({float(row["interval_sort"]) for row in rows})
    matrix = np.full((len(heights), len(intervals)), np.nan)
    for row in rows:
        h_index = heights.index(float(row["barrier_height"]))
        i_index = intervals.index(float(row["interval_sort"]))
        matrix[h_index, i_index] = float(row["unconditional_transmission_probability"])
    fig, ax = plt.subplots(figsize=(7, 4.5))
    image = ax.imshow(matrix, origin="lower", cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(intervals)), [f"{value:g}" for value in intervals])
    ax.set_yticks(np.arange(len(heights)), [f"{value:g}" for value in heights])
    ax.set_xlabel("measurement interval (0 = unmeasured)")
    ax.set_ylabel("barrier height")
    ax.set_title("No-click transmission heatmap")
    fig.colorbar(image, ax=ax, label="probability")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def _write_run_report(data: dict[str, Any], outputs: dict[str, Path], path: Path) -> Path:
    config = data["metadata"]["config"]
    lines = [
        f"# {config['name']} report",
        "",
        f"- Experiment: `{config['experiment']}`",
        f"- Grid: {config['grid']['n']} x {config['grid']['n']}, L={config['grid']['length']}",
        f"- Boundary: `{config['boundary']['kind']}`",
        f"- Measurement: `{config['measurement']['kind']}`",
        "",
        "## Figures",
        "",
        f"![Summary]({outputs['summary'].name})",
        "",
        f"![Potential]({outputs['potential'].name})",
        "",
    ]
    if "detected_probability" in data:
        lines.extend(
            [
                "## Final probabilities",
                "",
                "| Quantity | Value |",
                "| --- | ---: |",
                f"| Detector clicks | {_last(data.get('detected_probability')):.6f} |",
                f"| Absorbed at boundary | {_last(data.get('absorbed_probability')):.6f} |",
                f"| No-click survival | {_last(data.get('survival_weight')):.6f} |",
                f"| Conditional transmission | {_last(data.get('transmission_probability')):.6f} |",
                f"| Unconditional transmission | {_last(data.get('unconditional_transmission_probability')):.6f} |",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_comparison_report(
    rows: list[dict[str, float | str | None]], comparison: Path, heatmap: Path, path: Path
) -> Path:
    lines = [
        "# Zeno barrier sweep report",
        "",
        "This report compares no-click branch transmission, detector-click probability, boundary absorption, and no-click survival across barrier heights and repeated-measurement intervals.",
        "",
        "The measurement model treats each detector check as a no-click projection: probability in the detector region is accumulated as detector-click probability, removed from the conditional wavefunction, and the remaining branch is renormalized while survival weight records the unconditional no-click probability.",
        "",
        "Absorbing boundaries damp wavefunction amplitude near the edges after each split-step update; absorbed probability is tracked separately from detector clicks.",
        "",
        "The primary Zeno signal in this model is the unconditional transmission remaining in the no-click branch. Cumulative detector clicks are still reported, but they also reflect how many detector checks had a chance to remove already-transmitted amplitude.",
        "",
        "## Figures",
        "",
        f"![Comparison]({comparison.name})",
        "",
        f"![No-click transmission heatmap]({heatmap.name})",
        "",
        "## Sweep table",
        "",
        "| Run | V0 | interval | detector | absorbed | survival | unconditional transmission |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        interval = "unmeasured" if row["interval"] is None else f"{float(row['interval']):.6g}"
        lines.append(
            f"| {row['name']} | {float(row['barrier_height']):.6g} | {interval} | "
            f"{float(row['detected_probability']):.6f} | "
            f"{float(row['absorbed_probability']):.6f} | "
            f"{float(row['survival_weight']):.6f} | "
            f"{float(row['unconditional_transmission_probability']):.6f} |"
        )
    lines.extend(
        [
            "",
            "## Reading the sweep",
            "",
            "- Smaller measurement intervals correspond to more frequent detector checks.",
            "- A Zeno-like suppression signal appears when no-click branch transmission falls as detector checks become more frequent.",
            "- Absorbed probability is a boundary artifact control, not a detector event; large absorption indicates the simulation box or final time should be revisited.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
