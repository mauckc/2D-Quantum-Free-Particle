from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
import json
import tomllib

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .backends import create_backend
from .config import (
    BoundaryConfig,
    ExperimentConfig,
    GridConfig,
    MeasurementConfig,
    PotentialConfig,
    SolverConfig,
    WavePacketConfig,
)
from .experiments import run_experiment
from .grid import make_grid, normalize
from .solver import SplitStepSolver


@dataclass(frozen=True)
class NormCase:
    name: str
    n: int
    dt: float
    tf: float
    tolerance: float = 1e-8


@dataclass(frozen=True)
class DispersionCase:
    name: str
    n: int
    dt: float
    tf: float
    sigma: float
    tolerance: float = 0.35


@dataclass(frozen=True)
class BarrierCase:
    name: str
    heights: tuple[float, ...]
    tolerance: float = 1e-6


@dataclass(frozen=True)
class ZenoCase:
    name: str
    intervals: tuple[float, ...]
    tolerance: float = 1e-6


@dataclass(frozen=True)
class BackendCase:
    name: str
    backends: tuple[str, ...] = ("numpy", "pyfftw", "cupy")
    tolerance: float = 1e-10


@dataclass(frozen=True)
class ValidationConfig:
    name: str = "validation_suite"
    norm_cases: tuple[NormCase, ...] = field(default_factory=tuple)
    dispersion_cases: tuple[DispersionCase, ...] = field(default_factory=tuple)
    barrier_case: BarrierCase | None = None
    zeno_case: ZenoCase | None = None
    backend_case: BackendCase | None = None


def load_validation_config(path: str | Path) -> ValidationConfig:
    path = Path(path)
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    norm_cases = tuple(
        NormCase(
            name=str(item["name"]),
            n=int(item["n"]),
            dt=float(item["dt"]),
            tf=float(item["tf"]),
            tolerance=float(item.get("tolerance", NormCase.tolerance)),
        )
        for item in data.get("norm_cases", [])
    )
    dispersion_cases = tuple(
        DispersionCase(
            name=str(item["name"]),
            n=int(item["n"]),
            dt=float(item["dt"]),
            tf=float(item["tf"]),
            sigma=float(item["sigma"]),
            tolerance=float(item.get("tolerance", DispersionCase.tolerance)),
        )
        for item in data.get("dispersion_cases", [])
    )
    barrier_data = data.get("barrier_case")
    barrier_case = (
        None
        if barrier_data is None
        else BarrierCase(
            name=str(barrier_data.get("name", "barrier_height_sweep")),
            heights=tuple(float(value) for value in barrier_data["heights"]),
            tolerance=float(barrier_data.get("tolerance", BarrierCase.tolerance)),
        )
    )
    zeno_data = data.get("zeno_case")
    zeno_case = (
        None
        if zeno_data is None
        else ZenoCase(
            name=str(zeno_data.get("name", "zeno_interval_sweep")),
            intervals=tuple(float(value) for value in zeno_data["intervals"]),
            tolerance=float(zeno_data.get("tolerance", ZenoCase.tolerance)),
        )
    )
    backend_data = data.get("backend_case")
    backend_case = (
        None
        if backend_data is None
        else BackendCase(
            name=str(backend_data.get("name", "backend_parity")),
            backends=tuple(str(value) for value in backend_data.get("backends", BackendCase.backends)),
            tolerance=float(backend_data.get("tolerance", BackendCase.tolerance)),
        )
    )
    return ValidationConfig(
        name=str(data.get("name", path.stem)),
        norm_cases=norm_cases,
        dispersion_cases=dispersion_cases,
        barrier_case=barrier_case,
        zeno_case=zeno_case,
        backend_case=backend_case,
    )


def run_validation(config: ValidationConfig, out_dir: str | Path) -> dict[str, Any]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []

    for case in config.norm_cases:
        results.append(_validate_norm(case))
    for case in config.dispersion_cases:
        results.append(_validate_dispersion(case))
    if config.barrier_case is not None:
        results.append(_validate_barrier(config.barrier_case))
    if config.zeno_case is not None:
        results.append(_validate_zeno(config.zeno_case))
    if config.backend_case is not None:
        results.append(_validate_backends(config.backend_case))

    metrics = {
        "name": config.name,
        "passed": all(result["passed"] for result in results),
        "results": results,
        "config": asdict(config),
    }
    (out_dir / "validation_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8"
    )
    _plot_validation(metrics, out_dir / "validation_summary.png")
    _plot_convergence(metrics, out_dir)
    report = _write_validation_report(metrics, out_dir / "validation_report.md")
    _write_validation_html(report, out_dir / "validation_report.html")
    return metrics


def _validate_norm(case: NormCase) -> dict[str, Any]:
    result = run_experiment(
        ExperimentConfig(
            name=case.name,
            experiment="free_packet",
            grid=GridConfig(n=case.n, length=8.0),
            wave_packet=WavePacketConfig(center=(-1.5, 0.0), sigma=0.5, momentum=(2.0, 0.5)),
            potential=PotentialConfig(kind="free"),
            solver=SolverConfig(dt=case.dt, tf=case.tf, frame_interval=case.tf),
            boundary=BoundaryConfig(kind="periodic"),
        )
    )
    drift = float(np.max(np.abs(result.arrays["norm"] - 1.0)))
    return {
        "name": case.name,
        "kind": "norm_conservation",
        "value": drift,
        "tolerance": case.tolerance,
        "passed": drift <= case.tolerance,
        "details": {"n": case.n, "dt": case.dt, "tf": case.tf},
    }


def _validate_dispersion(case: DispersionCase) -> dict[str, Any]:
    result = run_experiment(
        ExperimentConfig(
            name=case.name,
            experiment="free_packet",
            grid=GridConfig(n=case.n, length=10.0),
            wave_packet=WavePacketConfig(center=(0.0, 0.0), sigma=case.sigma, momentum=(0.0, 0.0)),
            potential=PotentialConfig(kind="free"),
            solver=SolverConfig(dt=case.dt, tf=case.tf, frame_interval=case.tf),
            boundary=BoundaryConfig(kind="periodic"),
        )
    )
    probability = np.abs(result.arrays["psi_frames"][-1]) ** 2
    grid = result.grid
    variance_x = float(np.sum(probability * grid.X * grid.X) * grid.dx * grid.dy)
    initial_variance = case.sigma * case.sigma
    expected_variance = initial_variance + (case.tf * case.tf) / (4.0 * initial_variance)
    relative_error = abs(variance_x - expected_variance) / expected_variance
    return {
        "name": case.name,
        "kind": "free_gaussian_dispersion",
        "value": relative_error,
        "tolerance": case.tolerance,
        "passed": relative_error <= case.tolerance,
        "details": {
            "n": case.n,
            "dt": case.dt,
            "tf": case.tf,
            "sigma": case.sigma,
            "variance_x": variance_x,
            "expected_variance": expected_variance,
        },
    }


def _validate_barrier(case: BarrierCase) -> dict[str, Any]:
    transmissions = []
    for height in case.heights:
        result = run_experiment(_barrier_config(name=f"{case.name}_h{height:g}", height=height))
        transmissions.append(float(result.metrics["final_unconditional_transmission_probability"]))
    monotone = all(
        earlier + case.tolerance >= later
        for earlier, later in zip(transmissions, transmissions[1:])
    )
    return {
        "name": case.name,
        "kind": "barrier_transmission_monotonicity",
        "value": max(
            [0.0]
            + [later - earlier for earlier, later in zip(transmissions, transmissions[1:])]
        ),
        "tolerance": case.tolerance,
        "passed": monotone,
        "details": {"heights": list(case.heights), "transmissions": transmissions},
    }


def _validate_zeno(case: ZenoCase) -> dict[str, Any]:
    transmissions = []
    for interval in case.intervals:
        measurement = (
            MeasurementConfig(kind="none")
            if interval <= 0
            else MeasurementConfig(kind="zeno", interval=interval, detector_buffer=0.25)
        )
        result = run_experiment(
            _barrier_config(name=f"{case.name}_dt{interval:g}", height=40.0, measurement=measurement)
        )
        transmissions.append(float(result.metrics["final_unconditional_transmission_probability"]))
    unmeasured = transmissions[0]
    measured = transmissions[1:]
    suppresses = all(value <= unmeasured + case.tolerance for value in measured)
    frequent_first = all(
        earlier <= later + case.tolerance
        for earlier, later in zip(measured, measured[1:])
    )
    return {
        "name": case.name,
        "kind": "zeno_no_click_transmission",
        "value": max([0.0] + [value - unmeasured for value in measured]),
        "tolerance": case.tolerance,
        "passed": suppresses and frequent_first,
        "details": {"intervals": list(case.intervals), "transmissions": transmissions},
    }


def _validate_backends(case: BackendCase) -> dict[str, Any]:
    grid = make_grid(GridConfig(n=24, length=8.0))
    rng = np.random.default_rng(7)
    psi0 = normalize(
        (rng.normal(size=(grid.n, grid.n)) + 1j * rng.normal(size=(grid.n, grid.n))).astype(
            np.complex128
        ),
        grid,
    )
    potential = np.zeros((grid.n, grid.n), dtype=np.float64)
    solver_config = SolverConfig(dt=0.004, tf=0.004, frame_interval=0.004, backend="numpy")
    reference = SplitStepSolver(grid, potential, solver_config).step(psi0)
    comparisons = []
    for backend_name in case.backends:
        if backend_name == "numpy":
            continue
        try:
            backend = create_backend(backend_name)
        except RuntimeError as exc:
            comparisons.append({"backend": backend_name, "skipped": True, "reason": str(exc)})
            continue
        candidate = SplitStepSolver(grid, potential, solver_config, backend=backend).step(psi0)
        comparisons.append(
            {
                "backend": backend_name,
                "skipped": False,
                "max_abs_error": float(np.max(np.abs(candidate - reference))),
            }
        )
    active_errors = [
        item["max_abs_error"]
        for item in comparisons
        if not item.get("skipped", False)
    ]
    max_error = max([0.0] + active_errors)
    return {
        "name": case.name,
        "kind": "backend_parity",
        "value": max_error,
        "tolerance": case.tolerance,
        "passed": max_error <= case.tolerance,
        "details": {"comparisons": comparisons},
    }


def _barrier_config(
    name: str, height: float, measurement: MeasurementConfig | None = None
) -> ExperimentConfig:
    return ExperimentConfig(
        name=name,
        experiment="barrier_zeno",
        grid=GridConfig(n=40, length=8.0),
        wave_packet=WavePacketConfig(center=(-2.4, 0.0), sigma=0.45, momentum=(4.5, 0.0)),
        potential=PotentialConfig(kind="barrier", height=height, barrier_x=0.0, barrier_width=0.25),
        solver=SolverConfig(dt=0.002, tf=0.9, frame_interval=0.09),
        boundary=BoundaryConfig(kind="absorbing", width=1.0, strength=4.0),
        measurement=measurement or MeasurementConfig(kind="none"),
    )


def _plot_validation(metrics: dict[str, Any], path: Path) -> Path:
    results = metrics["results"]
    labels = [result["name"] for result in results]
    values = [float(result["value"]) for result in results]
    tolerances = [float(result["tolerance"]) for result in results]
    colors = ["#2ca25f" if result["passed"] else "#de2d26" for result in results]
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.2), 4.5))
    x = np.arange(len(labels))
    ax.bar(x, values, color=colors, label="observed")
    ax.scatter(x, tolerances, color="#1f2933", marker="_", s=140, label="tolerance")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylabel("validation metric")
    ax.set_title("Physics validation summary")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def _plot_convergence(metrics: dict[str, Any], out_dir: Path) -> None:
    norm_results = [
        result for result in metrics["results"] if result["kind"] == "norm_conservation"
    ]
    if norm_results:
        labels = [f"N={result['details']['n']}, dt={result['details']['dt']}" for result in norm_results]
        _plot_metric_series(
            labels,
            [float(result["value"]) for result in norm_results],
            [float(result["tolerance"]) for result in norm_results],
            "Norm drift across grid/time-step choices",
            "max |norm - 1|",
            out_dir / "norm_convergence.png",
        )

    dispersion_results = [
        result for result in metrics["results"] if result["kind"] == "free_gaussian_dispersion"
    ]
    if dispersion_results:
        labels = [
            f"N={result['details']['n']}, dt={result['details']['dt']}" for result in dispersion_results
        ]
        _plot_metric_series(
            labels,
            [float(result["value"]) for result in dispersion_results],
            [float(result["tolerance"]) for result in dispersion_results],
            "Free Gaussian dispersion error",
            "relative variance error",
            out_dir / "dispersion_error.png",
        )

    for result in metrics["results"]:
        if result["kind"] == "barrier_transmission_monotonicity":
            _plot_line_series(
                result["details"]["heights"],
                result["details"]["transmissions"],
                "Barrier transmission by potential height",
                "barrier height",
                "transmission probability",
                out_dir / "barrier_transmission.png",
            )
        if result["kind"] == "zeno_no_click_transmission":
            _plot_line_series(
                result["details"]["intervals"],
                result["details"]["transmissions"],
                "Zeno no-click transmission by measurement interval",
                "measurement interval",
                "unconditional transmission probability",
                out_dir / "zeno_transmission.png",
            )


def _plot_metric_series(
    labels: list[str],
    values: list[float],
    tolerances: list[float],
    title: str,
    ylabel: str,
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.4), 4.0))
    x = np.arange(len(labels))
    ax.plot(x, values, marker="o", color="#2b6cb0", label="observed")
    ax.plot(x, tolerances, marker="_", color="#9f1239", linestyle="none", markersize=12, label="tolerance")
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_line_series(
    x_values: list[float],
    y_values: list[float],
    title: str,
    xlabel: str,
    ylabel: str,
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.0))
    ax.plot(x_values, y_values, marker="o", color="#2f855a")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _write_validation_report(metrics: dict[str, Any], path: Path) -> Path:
    lines = [
        f"# {metrics['name']} validation report",
        "",
        f"Overall status: **{'PASS' if metrics['passed'] else 'FAIL'}**",
        "",
    ]
    for alt, filename in (
        ("Validation summary", "validation_summary.png"),
        ("Norm convergence", "norm_convergence.png"),
        ("Dispersion error", "dispersion_error.png"),
        ("Barrier transmission", "barrier_transmission.png"),
        ("Zeno transmission", "zeno_transmission.png"),
    ):
        if (path.parent / filename).exists():
            lines.extend([f"![{alt}]({filename})", ""])
    lines.extend(
        [
            "| Check | Kind | Value | Tolerance | Status |",
            "| --- | --- | ---: | ---: | --- |",
        ]
    )
    for result in metrics["results"]:
        status = "PASS" if result["passed"] else "FAIL"
        lines.append(
            f"| {result['name']} | {result['kind']} | {float(result['value']):.6g} | "
            f"{float(result['tolerance']):.6g} | {status} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- Norm checks report maximum deviation from unit total probability.",
            "- Dispersion checks compare measured Gaussian variance against the free-particle analytic variance.",
            "- Barrier checks require transmission to not increase as barrier height rises.",
            "- Zeno checks require measured no-click transmission to stay below the unmeasured baseline and increase as measurements become less frequent.",
            "- Backend checks skip unavailable optional packages and compare installed backends against NumPy.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_validation_html(markdown_path: Path, html_path: Path) -> Path:
    # Keep this local to avoid depending on the report renderer's private helpers.
    markdown = markdown_path.read_text(encoding="utf-8")
    lines = []
    for line in markdown.splitlines():
        if line.startswith("# "):
            lines.append(f"<h1>{line[2:]}</h1>")
        elif line.startswith("## "):
            lines.append(f"<h2>{line[3:]}</h2>")
        elif line.startswith("![") and "](" in line:
            alt, src = line[2:-1].split("](", 1)
            lines.append(f'<p><img src="{src}" alt="{alt}"></p>')
        elif line.startswith("|"):
            continue
        elif line.startswith("- "):
            lines.append(f"<li>{line[2:]}</li>")
        elif line.strip():
            lines.append(f"<p>{line}</p>")
    html_path.write_text(
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>Validation Report</title>"
        "<style>body{font-family:system-ui,sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;}"
        "img{max-width:100%;height:auto;border:1px solid #ddd}</style></head><body>"
        + "\n".join(lines)
        + "</body></html>",
        encoding="utf-8",
    )
    return html_path
