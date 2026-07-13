from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import json

import streamlit as st

from quantum_dynamics_lab.dashboard import (
    DashboardParameters,
    DashboardRun,
    DashboardSweep,
    experiment_summary,
    plot_diagnostic_panels,
    plot_double_slit_panels,
    plot_field_panels,
    run_dashboard_experiment,
    run_dashboard_validation,
    run_dashboard_zeno_sweep,
    timestamped_output_dir,
)

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    st.set_page_config(page_title="2D Quantum Dynamics Lab", layout="wide")
    st.title("2D Quantum Dynamics Lab")
    experiment_tab, sweep_tab, validation_tab = st.tabs(
        ["Experiment", "Zeno Sweep", "Validation"]
    )

    parameters, run_requested = _experiment_controls()
    with experiment_tab:
        if run_requested:
            out_dir = timestamped_output_dir(ROOT / "runs" / "dashboard", parameters.name)
            try:
                with st.spinner("Evolving wavefunction and generating artifacts..."):
                    st.session_state["dashboard_run"] = run_dashboard_experiment(
                        parameters, out_dir
                    )
            except (OSError, RuntimeError, ValueError) as exc:
                st.error(str(exc))
        dashboard_run = st.session_state.get("dashboard_run")
        if dashboard_run is None:
            st.info("Set experiment parameters in the sidebar and run a simulation.")
        else:
            _show_experiment(dashboard_run)

    with sweep_tab:
        _zeno_sweep_workspace(parameters)

    with validation_tab:
        _validation_workspace()


def _experiment_controls() -> tuple[DashboardParameters, bool]:
    with st.sidebar:
        st.header("Experiment setup")
        mode = st.segmented_control(
            "Experiment",
            ["Free packet", "Barrier / Zeno", "Double slit"],
            default="Barrier / Zeno",
            width="stretch",
        )
        is_barrier = mode == "Barrier / Zeno"
        is_double_slit = mode == "Double slit"
        experiment = (
            "double_slit"
            if is_double_slit
            else "barrier_zeno"
            if is_barrier
            else "free_packet"
        )
        mode_defaults = (
            {
                "grid_n": 64,
                "dt": 0.002,
                "tf": 1.2,
                "center_x": -3.2,
                "sigma": 0.45,
                "momentum_x": 6.0,
            }
            if is_double_slit
            else {
                "grid_n": 48,
                "dt": 0.004,
                "tf": 0.8,
                "center_x": -2.5,
                "sigma": 0.5,
                "momentum_x": 4.5,
            }
        )

        st.subheader("Grid and solver")
        grid_n = st.select_slider(
            "Grid points per axis",
            [32, 48, 64, 96],
            value=mode_defaults["grid_n"],
            key=f"{experiment}_grid_n",
        )
        grid_length = st.number_input(
            "Domain length",
            min_value=4.0,
            max_value=24.0,
            value=10.0,
            step=1.0,
            key=f"{experiment}_grid_length",
        )
        dt = st.number_input(
            "Time step",
            min_value=0.0005,
            max_value=0.02,
            value=mode_defaults["dt"],
            step=0.0005,
            format="%.4f",
            key=f"{experiment}_dt",
        )
        tf = st.number_input(
            "Final time",
            min_value=0.04,
            max_value=3.0,
            value=mode_defaults["tf"],
            step=0.1,
            key=f"{experiment}_tf",
        )
        frame_interval = st.number_input(
            "Saved-frame interval",
            min_value=float(dt),
            max_value=float(tf),
            value=min(0.04, float(tf)),
            step=float(dt),
            format="%.4f",
            key=f"{experiment}_frame_interval",
        )

        st.subheader("Wave packet")
        center_columns = st.columns(2)
        center_x = center_columns[0].number_input(
            "Center x",
            min_value=-10.0,
            max_value=10.0,
            value=mode_defaults["center_x"],
            step=0.25,
            key=f"{experiment}_center_x",
        )
        center_y = center_columns[1].number_input(
            "Center y",
            min_value=-10.0,
            max_value=10.0,
            value=0.0,
            step=0.25,
            key=f"{experiment}_center_y",
        )
        sigma = st.number_input(
            "Packet sigma",
            min_value=0.15,
            max_value=2.0,
            value=mode_defaults["sigma"],
            step=0.05,
            key=f"{experiment}_sigma",
        )
        momentum_columns = st.columns(2)
        momentum_x = momentum_columns[0].number_input(
            "Momentum x",
            min_value=-10.0,
            max_value=10.0,
            value=mode_defaults["momentum_x"],
            step=0.5,
            key=f"{experiment}_momentum_x",
        )
        momentum_y = momentum_columns[1].number_input(
            "Momentum y",
            min_value=-10.0,
            max_value=10.0,
            value=0.0,
            step=0.5,
            key=f"{experiment}_momentum_y",
        )

        barrier_height = 40.0
        barrier_x = 0.0
        barrier_width = 0.25
        measurement_kind = "none"
        measurement_interval: float | None = None
        measurement_time: float | None = None
        detector_buffer = 0.25
        slit_width = 0.45
        slit_separation = 1.2
        if is_barrier:
            st.subheader("Barrier and measurement")
            barrier_height = st.number_input(
                "Barrier height", min_value=0.0, max_value=200.0, value=40.0, step=5.0
            )
            barrier_columns = st.columns(2)
            barrier_x = barrier_columns[0].number_input(
                "Barrier x", min_value=-4.0, max_value=4.0, value=0.0, step=0.25
            )
            barrier_width = barrier_columns[1].number_input(
                "Barrier width", min_value=0.05, max_value=2.0, value=0.25, step=0.05
            )
            measurement_label = st.segmented_control(
                "Measurement", ["Unmeasured", "Zeno no-click"], default="Zeno no-click",
                width="stretch",
            )
            measurement_kind = "zeno" if measurement_label == "Zeno no-click" else "none"
            measurement_interval = st.number_input(
                "Measurement interval",
                min_value=float(dt),
                max_value=float(tf),
                value=min(0.12, float(tf)),
                step=float(dt),
                format="%.4f",
                disabled=measurement_kind == "none",
            )
            detector_buffer = st.number_input(
                "Detector buffer past barrier",
                min_value=0.0,
                max_value=3.0,
                value=0.25,
                step=0.05,
            )
        elif is_double_slit:
            st.subheader("Slit geometry and measurement")
            barrier_height = st.number_input(
                "Wall height",
                min_value=0.0,
                max_value=300.0,
                value=120.0,
                step=10.0,
            )
            barrier_columns = st.columns(2)
            barrier_x = barrier_columns[0].number_input(
                "Wall x", min_value=-4.0, max_value=4.0, value=0.0, step=0.25
            )
            barrier_width = barrier_columns[1].number_input(
                "Wall width", min_value=0.05, max_value=2.0, value=0.25, step=0.05
            )
            slit_columns = st.columns(2)
            slit_width = slit_columns[0].number_input(
                "Slit width", min_value=0.1, max_value=2.0, value=0.45, step=0.05
            )
            slit_separation = slit_columns[1].number_input(
                "Slit separation",
                min_value=0.2,
                max_value=5.0,
                value=1.2,
                step=0.1,
            )
            automatic_projection = st.toggle("Automatic projection time", value=True)
            explicit_time = st.number_input(
                "Which-path projection time",
                min_value=float(dt),
                max_value=float(tf),
                value=min(float(tf) * 0.55, float(tf) - float(dt)),
                step=float(dt),
                format="%.4f",
                disabled=automatic_projection,
            )
            measurement_kind = "which_path"
            measurement_time = None if automatic_projection else float(explicit_time)

        st.subheader("Boundary")
        absorbing_boundary = st.toggle("Absorbing edge", value=True)
        boundary_columns = st.columns(2)
        boundary_width = boundary_columns[0].number_input(
            "Absorber width", min_value=0.1, max_value=4.0, value=1.0, step=0.1,
            disabled=not absorbing_boundary,
        )
        boundary_strength = boundary_columns[1].number_input(
            "Absorber strength", min_value=0.1, max_value=20.0, value=4.0, step=0.5,
            disabled=not absorbing_boundary,
        )
        run_requested = st.button("Run simulation", type="primary", width="stretch")

    return (
        DashboardParameters(
            name=experiment,
            experiment=experiment,
            grid_n=int(grid_n),
            grid_length=float(grid_length),
            center_x=float(center_x),
            center_y=float(center_y),
            sigma=float(sigma),
            momentum_x=float(momentum_x),
            momentum_y=float(momentum_y),
            dt=float(dt),
            tf=float(tf),
            frame_interval=float(frame_interval),
            barrier_height=float(barrier_height),
            barrier_x=float(barrier_x),
            barrier_width=float(barrier_width),
            slit_width=float(slit_width),
            slit_separation=float(slit_separation),
            measurement_kind=measurement_kind,
            measurement_interval=(
                float(measurement_interval) if measurement_kind == "zeno" else None
            ),
            measurement_time=measurement_time,
            detector_buffer=float(detector_buffer),
            absorbing_boundary=absorbing_boundary,
            boundary_width=float(boundary_width),
            boundary_strength=float(boundary_strength),
        ),
        run_requested,
    )


def _show_experiment(dashboard_run: DashboardRun) -> None:
    result = dashboard_run.result
    metrics = experiment_summary(result)
    primary_metrics = st.columns(2)
    if result.config.experiment == "double_slit":
        primary_metrics[0].metric(
            "Final coherent norm", f"{float(metrics['final_norm']):.8f}"
        )
        primary_metrics[1].metric(
            "Final which-path norm", f"{float(metrics['final_which_path_norm']):.8f}"
        )
    else:
        primary_metrics[0].metric("Final norm", f"{float(metrics['final_norm']):.8f}")
        primary_metrics[1].metric(
            "Max norm drift", f"{float(metrics['max_norm_drift']):.2e}"
        )
    probability_metrics = st.columns(2)
    if result.config.experiment == "double_slit":
        probability_metrics[0].metric(
            "Coherent contrast",
            f"{float(metrics['coherent_interference_contrast']):.5f}",
        )
        probability_metrics[1].metric(
            "Which-path contrast",
            f"{float(metrics['which_path_interference_contrast']):.5f}",
        )
    elif result.config.experiment == "barrier_zeno":
        probability_metrics[0].metric(
            "No-click transmission",
            f"{float(metrics['final_unconditional_transmission_probability']):.5f}",
        )
        probability_metrics[1].metric(
            "Detector clicks", f"{float(metrics['final_detected_probability']):.5f}"
        )
    else:
        probability_metrics[0].metric(
            "Survival", f"{float(metrics['final_survival_weight']):.5f}"
        )
        probability_metrics[1].metric(
            "Absorbed", f"{float(metrics['final_absorbed_probability']):.5f}"
        )

    times = result.arrays["times"]
    frame_index = st.slider(
        "Saved frame",
        min_value=0,
        max_value=len(times) - 1,
        value=len(times) - 1,
        format="frame %d",
    )
    st.caption(f"t = {float(times[frame_index]):.4f}")
    field_columns = st.columns(3)
    field_figures = plot_field_panels(result, frame_index)
    for column, figure in zip(field_columns, field_figures, strict=True):
        column.pyplot(figure, width="stretch")
        plt.close(figure)

    if result.config.experiment == "double_slit":
        comparison_columns = st.columns(2)
        comparison_figures = plot_double_slit_panels(result, frame_index)
        for column, figure in zip(
            comparison_columns, comparison_figures, strict=True
        ):
            column.pyplot(figure, width="stretch")
            plt.close(figure)

    diagnostic_columns = st.columns(2)
    diagnostic_figures = plot_diagnostic_panels(result)
    for column, figure in zip(diagnostic_columns, diagnostic_figures, strict=True):
        column.pyplot(figure, width="stretch")
        plt.close(figure)

    st.subheader("Artifacts")
    st.code(str(dashboard_run.run_path.relative_to(ROOT)))
    download_columns = st.columns(3)
    download_columns[0].download_button(
        "Download run.npz",
        dashboard_run.run_path.read_bytes(),
        file_name="run.npz",
        mime="application/octet-stream",
        width="stretch",
        on_click="ignore",
    )
    metrics_path = dashboard_run.run_path.parent / "metrics.json"
    download_columns[1].download_button(
        "Download metrics.json",
        metrics_path.read_bytes(),
        file_name="metrics.json",
        mime="application/json",
        width="stretch",
        on_click="ignore",
    )
    report_path = dashboard_run.outputs["html_report"]
    download_columns[2].download_button(
        "Download report.html",
        report_path.read_bytes(),
        file_name="report.html",
        mime="text/html",
        width="stretch",
        on_click="ignore",
    )


def _zeno_sweep_workspace(parameters: DashboardParameters) -> None:
    interval_options = [
        value
        for value in [0.0, 0.04, 0.06, 0.08, 0.12, 0.24]
        if value == 0 or value <= parameters.tf
    ]
    default_intervals = [
        value for value in [0.0, 0.06, 0.12] if value in interval_options
    ]
    if len(default_intervals) == 1 and len(interval_options) > 1:
        default_intervals.append(interval_options[1])
    control_columns = st.columns(2)
    heights = control_columns[0].multiselect(
        "Barrier heights",
        options=[20.0, 40.0, 60.0, 80.0, 120.0],
        default=[20.0, 40.0],
        max_selections=5,
        format_func=lambda value: f"{value:g}",
    )
    intervals = control_columns[1].multiselect(
        "Measurement intervals",
        options=interval_options,
        default=default_intervals,
        max_selections=5,
        format_func=lambda value: "Unmeasured" if value == 0 else f"{value:g}",
    )
    variant_count = len(heights) * len(intervals)
    st.caption(f"{variant_count} variants using the sidebar grid, packet, and solver settings")
    if st.button(
        "Run Zeno sweep",
        type="primary",
        disabled=(
            variant_count == 0
            or variant_count > 25
            or 0.0 not in intervals
            or not any(value > 0 for value in intervals)
        ),
    ):
        out_dir = timestamped_output_dir(
            ROOT / "reports" / "dashboard", "zeno-sweep"
        )
        base_parameters = replace(
            parameters,
            name="dashboard_zeno_sweep",
            experiment="barrier_zeno",
            measurement_kind="none",
            measurement_interval=None,
            measurement_time=None,
        )
        try:
            with st.spinner("Running controlled Zeno variants and building reports..."):
                st.session_state["dashboard_zeno_sweep"] = run_dashboard_zeno_sweep(
                    base_parameters,
                    tuple(float(value) for value in heights),
                    tuple(float(value) for value in intervals),
                    out_dir,
                )
        except (OSError, RuntimeError, ValueError) as exc:
            st.error(str(exc))

    sweep = st.session_state.get("dashboard_zeno_sweep")
    if sweep is None:
        st.info("Select a controlled parameter grid and run the sweep.")
        return
    _show_zeno_sweep(sweep)


def _show_zeno_sweep(sweep: DashboardSweep) -> None:
    plot_columns = st.columns(2)
    plot_columns[0].image(sweep.outputs["comparison"], width="stretch")
    plot_columns[1].image(sweep.outputs["heatmap"], width="stretch")
    st.dataframe(
        list(sweep.rows),
        width="stretch",
        hide_index=True,
        column_config={
            "barrier_height": st.column_config.NumberColumn(format="%.1f"),
            "measurement_interval": st.column_config.NumberColumn(format="%.3f"),
            "transmission": st.column_config.NumberColumn(format="%.3e"),
            "detected": st.column_config.NumberColumn(format="%.3e"),
            "survival": st.column_config.NumberColumn(format="%.3e"),
            "absorbed": st.column_config.NumberColumn(format="%.3e"),
        },
    )
    st.code(str(sweep.out_dir.relative_to(ROOT)))
    download_columns = st.columns(2)
    download_columns[0].download_button(
        "Download comparison.csv",
        sweep.outputs["csv"].read_bytes(),
        file_name="comparison.csv",
        mime="text/csv",
        width="stretch",
        on_click="ignore",
    )
    download_columns[1].download_button(
        "Download sweep report",
        sweep.outputs["html_report"].read_bytes(),
        file_name="report.html",
        mime="text/html",
        width="stretch",
        on_click="ignore",
    )


def _validation_workspace() -> None:
    config_value = st.text_input(
        "Validation config", value="examples/validation_suite.toml"
    )
    if st.button("Run validation suite", type="primary"):
        config_path = Path(config_value)
        if not config_path.is_absolute():
            config_path = ROOT / config_path
        out_dir = timestamped_output_dir(
            ROOT / "reports" / "dashboard", "validation-suite"
        )
        try:
            with st.spinner("Running deterministic physics benchmarks..."):
                metrics = run_dashboard_validation(config_path, out_dir)
            st.session_state["dashboard_validation"] = (metrics, out_dir)
        except (OSError, RuntimeError, ValueError) as exc:
            st.error(str(exc))

    validation_state = st.session_state.get("dashboard_validation")
    if validation_state is None:
        st.info("Run the validation suite to inspect benchmark status and convergence plots.")
        return
    metrics, out_dir = validation_state
    if metrics["passed"]:
        st.success("All configured validation checks passed.")
    else:
        st.error("One or more validation checks failed.")
    rows = [
        {
            "check": result["name"],
            "kind": result["kind"],
            "value": float(result["value"]),
            "tolerance": float(result["tolerance"]),
            "status": "PASS" if result["passed"] else "FAIL",
        }
        for result in metrics["results"]
    ]
    st.dataframe(
        rows,
        width="stretch",
        hide_index=True,
        column_config={
            "value": st.column_config.NumberColumn(format="%.3e"),
            "tolerance": st.column_config.NumberColumn(format="%.3e"),
        },
    )
    plot_names = (
        "validation_summary.png",
        "norm_convergence.png",
        "dispersion_error.png",
        "barrier_transmission.png",
        "zeno_transmission.png",
    )
    plots = [out_dir / name for name in plot_names if (out_dir / name).exists()]
    for start in range(0, len(plots), 2):
        plot_columns = st.columns(2)
        for column, path in zip(plot_columns, plots[start : start + 2]):
            column.image(
                path,
                caption=path.stem.replace("_", " "),
                width="stretch",
            )
    st.code(str(out_dir.relative_to(ROOT)))
    validation_json = out_dir / "validation_metrics.json"
    st.download_button(
        "Download validation metrics",
        validation_json.read_bytes(),
        file_name="validation_metrics.json",
        mime="application/json",
        on_click="ignore",
    )
    with st.expander("Benchmark details"):
        st.json(json.loads(validation_json.read_text(encoding="utf-8")))


if __name__ == "__main__":
    main()
