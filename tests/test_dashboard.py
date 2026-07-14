from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from quantum_dynamics_lab.dashboard import (
    DashboardParameters,
    build_experiment_config,
    experiment_summary,
    load_flagship_dashboard_summary,
    plot_diagnostic_panels,
    plot_double_slit_panels,
    plot_field_panels,
    plot_wave_field_panels,
    run_dashboard_experiment,
    run_dashboard_wave_forward,
    run_dashboard_zeno_sweep,
)


def test_dashboard_builds_free_and_zeno_configs() -> None:
    free = build_experiment_config(
        DashboardParameters(
            name="free",
            experiment="free_packet",
            measurement_kind="zeno",
            measurement_interval=0.04,
        )
    )
    zeno = build_experiment_config(DashboardParameters())

    assert free.potential.kind == "free"
    assert free.measurement.kind == "none"
    assert zeno.potential.kind == "barrier"
    assert zeno.measurement.kind == "zeno"
    assert zeno.output.render_video is False


def test_dashboard_rejects_detector_inside_absorber() -> None:
    parameters = replace(
        DashboardParameters(), grid_length=8.0, boundary_width=2.0, detector_buffer=2.0
    )

    with pytest.raises(ValueError, match="detector must lie before"):
        build_experiment_config(parameters)


def test_dashboard_rejects_barrier_crossing_grid_edge() -> None:
    parameters = replace(
        DashboardParameters(), barrier_x=4.9, barrier_width=0.4, absorbing_boundary=False
    )

    with pytest.raises(ValueError, match="barrier must fit entirely"):
        build_experiment_config(parameters)


def test_dashboard_builds_double_slit_config() -> None:
    config = build_experiment_config(
        DashboardParameters(
            name="double",
            experiment="double_slit",
            measurement_kind="which_path",
            measurement_interval=None,
            measurement_time=0.4,
        )
    )

    assert config.potential.kind == "double_slit"
    assert config.potential.slit_width == 0.45
    assert config.measurement.kind == "which_path"
    assert config.measurement.time == 0.4


def test_dashboard_rejects_overlapping_slits() -> None:
    parameters = replace(
        DashboardParameters(),
        experiment="double_slit",
        measurement_kind="which_path",
        measurement_interval=None,
        slit_width=0.8,
        slit_separation=0.7,
    )

    with pytest.raises(ValueError, match="separation must exceed"):
        build_experiment_config(parameters)


def test_dashboard_run_saves_artifacts_and_plots(tmp_path) -> None:
    parameters = DashboardParameters(
        name="dashboard_smoke",
        experiment="barrier_zeno",
        grid_n=24,
        grid_length=8.0,
        center_x=-2.0,
        sigma=0.5,
        momentum_x=4.0,
        dt=0.004,
        tf=0.04,
        frame_interval=0.02,
        measurement_kind="zeno",
        measurement_interval=0.02,
        boundary_width=1.0,
    )

    dashboard_run = run_dashboard_experiment(parameters, tmp_path / "run")

    assert dashboard_run.run_path.exists()
    assert (dashboard_run.run_path.parent / "metrics.json").exists()
    assert dashboard_run.outputs["summary"].exists()
    assert dashboard_run.outputs["html_report"].exists()
    summary = experiment_summary(dashboard_run.result)
    assert summary["max_norm_drift"] < 1e-8
    field_figures = plot_field_panels(dashboard_run.result, 0)
    diagnostic_figures = plot_diagnostic_panels(dashboard_run.result)
    assert len(field_figures) == 3
    assert all(len(figure.axes) == 2 for figure in field_figures)
    assert len(diagnostic_figures) == 2
    assert all(len(figure.axes) == 1 for figure in diagnostic_figures)
    for figure in (*field_figures, *diagnostic_figures):
        plt.close(figure)


def test_dashboard_double_slit_run_builds_comparison_panels(tmp_path) -> None:
    parameters = DashboardParameters(
        name="double_slit_smoke",
        experiment="double_slit",
        grid_n=24,
        grid_length=8.0,
        center_x=-2.0,
        sigma=0.5,
        momentum_x=4.0,
        dt=0.004,
        tf=0.08,
        frame_interval=0.04,
        measurement_kind="which_path",
        measurement_interval=None,
        measurement_time=0.04,
        boundary_width=1.0,
    )

    dashboard_run = run_dashboard_experiment(parameters, tmp_path / "double")
    summary = experiment_summary(dashboard_run.result)
    comparison_figures = plot_double_slit_panels(dashboard_run.result, 0)
    diagnostic_figures = plot_diagnostic_panels(dashboard_run.result)

    assert dashboard_run.run_path.exists()
    assert "coherent_interference_contrast" in summary
    assert "which_path_interference_contrast" in summary
    assert len(comparison_figures) == 2
    assert len(diagnostic_figures) == 2
    for figure in (*comparison_figures, *diagnostic_figures):
        plt.close(figure)


def test_dashboard_zeno_sweep_writes_comparison_artifacts(tmp_path) -> None:
    parameters = DashboardParameters(
        name="sweep_smoke",
        experiment="barrier_zeno",
        grid_n=24,
        grid_length=8.0,
        center_x=-2.0,
        sigma=0.5,
        momentum_x=4.0,
        dt=0.004,
        tf=0.08,
        frame_interval=0.04,
        measurement_kind="none",
        measurement_interval=None,
        boundary_width=1.0,
    )

    sweep = run_dashboard_zeno_sweep(
        parameters, (20.0,), (0.0, 0.04), tmp_path / "sweep"
    )

    assert len(sweep.run_paths) == 2
    assert len(sweep.rows) == 2
    assert sweep.outputs["comparison"].exists()
    assert sweep.outputs["heatmap"].exists()
    assert sweep.outputs["csv"].exists()
    assert sweep.outputs["html_report"].exists()


def test_dashboard_wave_forward_and_flagship_reference(tmp_path) -> None:
    repository_root = Path(__file__).parents[1]
    wave = run_dashboard_wave_forward(
        repository_root / "examples" / "gaussian_lens.toml",
        tmp_path / "wave",
    )
    figures = plot_wave_field_panels(wave.run)
    reference = load_flagship_dashboard_summary(repository_root)

    assert wave.run_path.exists()
    assert wave.run.metrics["power_efficiency"] == pytest.approx(1.0, abs=1e-13)
    assert len(figures) == 2
    assert all(len(figure.axes) == 4 for figure in figures)
    assert reference.comparison["robustness_claim_passed"]
    assert reference.validation["trusted"]
    for figure in figures:
        plt.close(figure)


def test_streamlit_dashboard_starts_without_errors() -> None:
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    app_path = Path(__file__).parents[1] / "apps" / "streamlit_dashboard.py"
    app = AppTest.from_file(str(app_path), default_timeout=10).run()

    assert not app.exception
    assert app.title[0].value == "Differentiable 2D Wave Inverse-Design Lab"
    assert [tab.label for tab in app.tabs[:4]] == [
        "Forward Optics",
        "Inverse Design",
        "Validated Flagship",
        "Legacy Quantum",
    ]
