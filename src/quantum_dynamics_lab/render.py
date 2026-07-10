from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import html
import json
import shutil

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np


def render_research_narrative(
    zeno_dir: str | Path, double_slit_run: str | Path, out_dir: str | Path
) -> dict[str, Path]:
    from .artifacts import load_run

    zeno_dir = Path(zeno_dir)
    double_slit_run = Path(double_slit_run)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = _read_comparison_csv(zeno_dir / "comparison.csv")
    double_data = load_run(double_slit_run)
    double_metrics_path = double_slit_run.parent / "metrics.json"
    double_metrics = (
        json.loads(double_metrics_path.read_text(encoding="utf-8"))
        if double_metrics_path.exists()
        else {}
    )
    zeno_comparison = _copy_if_exists(zeno_dir / "comparison.png", out_dir / "zeno_comparison.png")
    zeno_heatmap = _copy_if_exists(
        zeno_dir / "zeno_transmission_heatmap.png",
        out_dir / "zeno_transmission_heatmap.png",
    )
    double_summary = _copy_if_exists(
        double_slit_run.parent / "summary.png", out_dir / "double_slit_summary.png"
    )

    narrative_path = _write_full_research_narrative(
        rows=rows,
        zeno_comparison=zeno_comparison,
        zeno_heatmap=zeno_heatmap,
        double_summary=double_summary,
        double_config=double_data["metadata"]["config"],
        double_metrics=double_metrics,
        path=out_dir / "research_narrative.md",
    )
    html_path = _write_html_from_markdown(narrative_path, out_dir / "research_narrative.html")
    return {"narrative": narrative_path, "html": html_path}


def _copy_if_exists(source: Path, destination: Path) -> Path:
    if source.exists():
        shutil.copy2(source, destination)
    return destination


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
    outputs["html_report"] = _write_html_from_markdown(
        outputs["report"], out_dir / "report.html"
    )
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
    report_path = _write_comparison_report(rows, path, heatmap_path, out_dir / "report.md")
    _write_html_from_markdown(report_path, out_dir / "report.html")
    _write_explorer(rows, out_dir / "explorer.html")
    _write_research_narrative(rows, path, heatmap_path, out_dir / "research_narrative.md")
    _write_html_from_markdown(out_dir / "research_narrative.md", out_dir / "research_narrative.html")
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


def _read_comparison_csv(path: Path) -> list[dict[str, float | str | None]]:
    rows: list[dict[str, float | str | None]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, float | str | None] = {"name": row["name"]}
            for key, value in row.items():
                if key == "name":
                    continue
                if value == "":
                    parsed[key] = None
                else:
                    parsed[key] = float(value)
            rows.append(parsed)
    return rows


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


def _write_html_from_markdown(markdown_path: Path, html_path: Path) -> Path:
    markdown = markdown_path.read_text(encoding="utf-8")
    body = _simple_markdown_to_html(markdown)
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(markdown_path.stem.replace('_', ' ').title())}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.5; max-width: 1080px; margin: 2rem auto; padding: 0 1rem; color: #17202a; }}
    h1, h2 {{ line-height: 1.2; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #d6dbe0; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #d6dbe0; padding: 0.45rem 0.6rem; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
    code {{ background: #f1f3f5; padding: 0.1rem 0.25rem; border-radius: 4px; }}
  </style>
</head>
<body>
{body}
</body>
</html>
"""
    html_path.write_text(html_text, encoding="utf-8")
    return html_path


def _simple_markdown_to_html(markdown: str) -> str:
    lines = markdown.splitlines()
    output: list[str] = []
    table: list[str] = []
    in_list = False
    for line in lines:
        if table and not line.startswith("|"):
            output.append(_table_to_html(table))
            table = []
        if in_list and not line.startswith("- "):
            output.append("</ul>")
            in_list = False
        if line.startswith("|"):
            table.append(line)
        elif line.startswith("# "):
            output.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            output.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("![") and "](" in line and line.endswith(")"):
            alt, src = line[2:-1].split("](", 1)
            output.append(f'<p><img src="{html.escape(src)}" alt="{html.escape(alt)}"></p>')
        elif line.startswith("- "):
            if not in_list:
                output.append("<ul>")
                in_list = True
            output.append(f"<li>{_inline_markdown(line[2:])}</li>")
        elif not line.strip():
            output.append("")
        else:
            output.append(f"<p>{_inline_markdown(line)}</p>")
    if table:
        output.append(_table_to_html(table))
    if in_list:
        output.append("</ul>")
    return "\n".join(output)


def _table_to_html(lines: list[str]) -> str:
    rows = []
    for index, line in enumerate(lines):
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if index == 1 and all(set(cell) <= {"-", ":"} for cell in cells):
            continue
        tag = "th" if index == 0 else "td"
        rows.append("<tr>" + "".join(f"<{tag}>{_inline_markdown(cell)}</{tag}>" for cell in cells) + "</tr>")
    return "<table>\n" + "\n".join(rows) + "\n</table>"


def _inline_markdown(text: str) -> str:
    escaped = html.escape(text)
    parts = escaped.split("`")
    for index in range(1, len(parts), 2):
        parts[index] = f"<code>{parts[index]}</code>"
    return "".join(parts)


def _write_explorer(rows: list[dict[str, float | str | None]], path: Path) -> Path:
    data = json.dumps(rows)
    heights = sorted({float(row["barrier_height"]) for row in rows})
    intervals = sorted({float(row["interval_sort"]) for row in rows})
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Zeno Parameter Explorer</title>
  <style>
    body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem auto; max-width: 980px; padding: 0 1rem; color: #17202a; }}
    .controls {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; margin: 1.5rem 0; }}
    label {{ display: grid; gap: 0.35rem; font-weight: 600; }}
    select {{ font: inherit; padding: 0.4rem; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 0.75rem; }}
    .metric {{ border: 1px solid #d6dbe0; padding: 0.8rem; }}
    .metric span {{ display: block; font-size: 0.85rem; color: #52616f; }}
    .metric strong {{ font-size: 1.35rem; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 1.5rem; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #d6dbe0; padding: 0.45rem; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
  </style>
</head>
<body>
  <h1>Zeno Parameter Explorer</h1>
  <p>Inspect how barrier height and measurement interval change no-click transmission, detector clicks, survival weight, and boundary absorption.</p>
  <div class="controls">
    <label>Barrier height<select id="height"></select></label>
    <label>Measurement interval<select id="interval"></select></label>
  </div>
  <div class="metrics" id="metrics"></div>
  <table id="table"></table>
  <script>
    const rows = {data};
    const heights = {json.dumps(heights)};
    const intervals = {json.dumps(intervals)};
    const heightSelect = document.querySelector("#height");
    const intervalSelect = document.querySelector("#interval");
    const metrics = document.querySelector("#metrics");
    const table = document.querySelector("#table");
    for (const value of heights) heightSelect.add(new Option(value, value));
    for (const value of intervals) intervalSelect.add(new Option(value === 0 ? "unmeasured" : value, value));
    function selectedRow() {{
      return rows.find(row => Number(row.barrier_height) === Number(heightSelect.value) && Number(row.interval_sort) === Number(intervalSelect.value));
    }}
    function fmt(value) {{ return Number(value).toFixed(6); }}
    function render() {{
      const row = selectedRow();
      metrics.innerHTML = "";
      for (const [label, key] of [
        ["No-click transmission", "unconditional_transmission_probability"],
        ["Detector clicks", "detected_probability"],
        ["No-click survival", "survival_weight"],
        ["Boundary absorption", "absorbed_probability"],
      ]) {{
        const div = document.createElement("div");
        div.className = "metric";
        div.innerHTML = `<span>${{label}}</span><strong>${{fmt(row[key])}}</strong>`;
        metrics.appendChild(div);
      }}
      table.innerHTML = `<tr><th>Run</th><th>V0</th><th>interval</th><th>detector x</th><th>unconditional transmission</th><th>detector clicks</th><th>absorbed</th><th>survival</th></tr>` +
        rows.map(r => `<tr><td>${{r.name}}</td><td>${{r.barrier_height}}</td><td>${{r.interval ?? "unmeasured"}}</td><td>${{fmt(r.detector_x)}}</td><td>${{fmt(r.unconditional_transmission_probability)}}</td><td>${{fmt(r.detected_probability)}}</td><td>${{fmt(r.absorbed_probability)}}</td><td>${{fmt(r.survival_weight)}}</td></tr>`).join("");
    }}
    heightSelect.addEventListener("change", render);
    intervalSelect.addEventListener("change", render);
    render();
  </script>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")
    return path


def _write_research_narrative(
    rows: list[dict[str, float | str | None]], comparison: Path, heatmap: Path, path: Path
) -> Path:
    unmeasured = [row for row in rows if row["interval"] is None]
    measured = [row for row in rows if row["interval"] is not None]
    strongest = min(measured, key=lambda row: float(row["unconditional_transmission_probability"]))
    weakest_barrier = max(unmeasured, key=lambda row: float(row["unconditional_transmission_probability"]))
    lines = [
        "# Reproducible Zeno Barrier Narrative",
        "",
        "## Question",
        "",
        "Can repeated no-click observation suppress transmission of a Gaussian wave packet through a finite barrier on a two-dimensional lattice?",
        "",
        "## Method",
        "",
        "Each run evolves the same packet with a split-step Fourier method. The sweep varies barrier height and measurement interval. Measurements are modeled as no-click projections at a detector region to the right of the barrier, while absorbing boundaries track probability that reaches the edge of the simulation box.",
        "",
        "## Figures",
        "",
        f"![Comparison]({comparison.name})",
        "",
        f"![No-click transmission heatmap]({heatmap.name})",
        "",
        "## Result",
        "",
        f"The largest unmeasured no-click transmission in this sweep is {float(weakest_barrier['unconditional_transmission_probability']):.6f} for `{weakest_barrier['name']}`. The strongest suppression among measured runs leaves {float(strongest['unconditional_transmission_probability']):.6f} unconditional no-click transmission for `{strongest['name']}`.",
        "",
        "## Assumptions",
        "",
        "- Periodic FFT evolution is combined with an optional absorbing edge mask to reduce wraparound artifacts.",
        "- Detector clicks are accumulated separately from boundary absorption.",
        "- The Zeno signal is interpreted through unconditional no-click branch transmission, not cumulative detector-click probability alone.",
        "- The generated figures and CSV are reproducible from the committed TOML configs and the `quantum-lab compare` command.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_full_research_narrative(
    rows: list[dict[str, float | str | None]],
    zeno_comparison: Path,
    zeno_heatmap: Path,
    double_summary: Path,
    double_config: dict[str, Any],
    double_metrics: dict[str, Any],
    path: Path,
) -> Path:
    measured = [row for row in rows if row["interval"] is not None]
    unmeasured = [row for row in rows if row["interval"] is None]
    strongest = min(measured, key=lambda row: float(row["unconditional_transmission_probability"]))
    baseline = max(unmeasured, key=lambda row: float(row["unconditional_transmission_probability"]))
    coherent = float(double_metrics.get("coherent_interference_contrast", 0.0))
    which_path = float(double_metrics.get("which_path_interference_contrast", 0.0))
    lines = [
        "# 2D Quantum Dynamics Research Narrative",
        "",
        "## Research Questions",
        "",
        "1. Does repeated no-click observation suppress barrier transmission on a two-dimensional lattice?",
        "2. Does which-path projection reduce double-slit interference contrast?",
        "",
        "## Reproducible Inputs",
        "",
        "- Zeno sweep CSV: `comparison.csv`",
        f"- Double-slit run: `{double_config['name']}`",
        f"- Double-slit grid: {double_config['grid']['n']} x {double_config['grid']['n']}, L={double_config['grid']['length']}",
        "",
        "## Zeno Barrier Result",
        "",
        f"The unmeasured baseline with the largest transmission is `{baseline['name']}` with unconditional transmission {float(baseline['unconditional_transmission_probability']):.6f}. The strongest measured suppression is `{strongest['name']}` with unconditional no-click transmission {float(strongest['unconditional_transmission_probability']):.6f}.",
        "",
        f"![Zeno comparison]({zeno_comparison.name})",
        "",
        f"![Zeno heatmap]({zeno_heatmap.name})",
        "",
        "## Double-Slit Result",
        "",
        f"The coherent screen-profile contrast is {coherent:.6f}; the which-path contrast is {which_path:.6f}. Lower which-path contrast is consistent with measurement suppressing interference structure in the branch-combined probability.",
        "",
        f"![Double-slit summary]({double_summary.name})",
        "",
        "## Assumptions And Limits",
        "",
        "- The solver uses a split-step Fourier method and optional FFT backends; numerical behavior should be checked when changing backend or grid resolution.",
        "- Absorbing edges reduce wraparound but introduce a boundary-loss diagnostic that should remain small enough for the intended interpretation.",
        "- The Zeno analysis focuses on unconditional no-click branch transmission rather than raw cumulative detector clicks.",
        "- The double-slit comparison is an operational which-path projection model, not a full detector-environment model.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
