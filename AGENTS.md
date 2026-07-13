# Project Operating Context

## Mission

This repository is being incrementally pivoted from a Quantum Zeno demonstration
into a differentiable 2D wave inverse-design lab, with scalar photonics as the
first practical application. The existing quantum workflows remain supported as
legacy demonstrations and regression tests while a generic propagation core is
built underneath them.

## Read Before Making Changes

1. Read `ROADMAP.md` for the current milestone, scope, and exit criteria.
2. Read the accepted decisions in `docs/adr/` that apply to the work.
3. Follow `CONTRIBUTING.md` for issue, branch, PR, validation, and evidence rules.

GitHub issues are the source of truth for active work. `ROADMAP.md` records stable
direction and milestone gates; it is not a daily task list.

## Current Priority

**M0 - Pivot foundation** is complete. The next behavioral work is the first
part of **M1 - Forward optics v0.2**: extract a pure, backend-neutral
propagation kernel and then make the maintained quantum solver delegate to it
without changing recorded quantum v0.1 results. Do not add JAX, optimization,
or more dashboard features ahead of that compatibility-preserving refactor.

## Project Invariants

- Do not remove or silently break the existing `quantum-lab` CLI, committed
  quantum configs, or dashboard workflows during the pivot.
- Do not present the current no-click Zeno trend or global min/max double-slit
  contrast as independently validated research evidence.
- Keep numerical propagation separate from configuration parsing, plotting,
  artifact writing, dashboards, and experiment-specific branching.
- Keep NumPy as the correctness reference. Add JAX as the differentiable path,
  not as an implicit replacement for the reference implementation.
- The first optics release is monochromatic, coherent, scalar free-space optics.
  Do not imply vector-Maxwell or high-index integrated-photonics accuracy.
- Generated `runs/` and `reports/` remain untracked. Commit inputs, small
  reference metrics, and regeneration code instead of large outputs.

## Required Evidence

Match evidence to the change:

- Refactors: legacy regression and backend-parity tests.
- Propagation models: analytic checks, conservation/reversibility, and grid
  convergence.
- Differentiation: centered directional finite-difference gradient checks.
- Optimization: deterministic objective-improvement tests plus saved history.
- Research results: declared baselines, held-out perturbations, and a generated
  model-transfer report.

Run at least `uv run pytest` before handing off code changes. Add narrower
commands and generated evidence when required by the issue acceptance criteria.
