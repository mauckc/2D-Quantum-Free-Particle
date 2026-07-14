# Differentiable 2D Wave Inverse-Design Lab Roadmap

Last updated: 2026-07-13

## Direction

The project will evolve from a maintained 2D quantum-dynamics demonstration into
a reproducible wave-propagation and inverse-design laboratory. Scalar photonics
is the first application because its FFT propagation operators closely match the
existing split-step solver.

The pivot is incremental. Existing quantum examples stay runnable while their
numerical core is generalized and reused by a new optics layer.

## First Product Boundary

The first inverse-design release will support:

- monochromatic coherent scalar fields `E(x, y; z)`;
- NumPy reference propagation;
- JAX-differentiable propagation;
- Gaussian, Hermite-Gaussian, and arbitrary complex input fields;
- Fresnel/paraxial and band-limited angular-spectrum propagation;
- apertures, thin lenses, and phase-only masks;
- mode-overlap, intensity-matching, and power-efficiency objectives;
- bounded phase, smoothing, quantization, and robustness constraints;
- reproducible configs, metrics, optimization histories, plots, and reports.

The first release will not claim support for:

- polarization or vector Maxwell fields;
- reliable high-index-contrast integrated photonics;
- nonlinear, broadband, or partially coherent propagation;
- manufacturing-ready layout or mask export.

## Delivery Flow

```text
M0 Pivot foundation
  -> M1 Forward optics v0.2
  -> M2 Differentiable solver v0.3
  -> M3 Inverse-design MVP v0.4
  -> M4 Robust flagship v1.0
```

## Milestones

| Milestone | Status | Outcome |
| --- | --- | --- |
| M0 - Pivot foundation | Complete | Durable context, workflow templates, accepted decisions, corrected framing, and a frozen quantum v0.1 regression baseline |
| M1 - Forward optics v0.2 | Complete | Pure propagation core plus validated scalar-optics sources, elements, Fresnel propagation, and angular-spectrum propagation |
| M2 - Differentiable solver v0.3 | Complete | JAX parity, JIT-compatible propagation, mode-overlap objectives, and validated gradients |
| M3 - Inverse-design MVP v0.4 | Planned | Constrained phase parameterization, regularization, optimization runner, checkpoints, and deterministic smoke optimization |
| M4 - Robust flagship v1.0 | Planned | Robust multi-plane mode conversion, held-out perturbation results, cross-model validation, and the new CLI/dashboard identity |

### M0 - Pivot Foundation

Deliverables:

- repository-owned roadmap and agent context;
- architecture decision records;
- issue and pull-request templates;
- contribution, validation, and evidence rules;
- corrected README language for the maintained algorithm and scientific claims;
- small quantum regression fixtures covering free propagation, barrier runs,
  and double-slit outputs;
- a tagged or otherwise recorded v0.1 baseline after the above is reviewed.

Exit gate:

- all current CLI and dashboard workflows remain available;
- the existing test suite passes;
- future work can be understood from a fresh checkout without this conversation;
- the first M1 issues have explicit acceptance criteria.

### M1 - Forward Optics v0.2

Deliverables:

- a pure propagation API with no plotting, I/O, or experiment branching;
- the existing quantum solver delegating to the new core;
- optics-specific grids, units, sources, targets, and elements;
- Fresnel/paraxial and band-limited angular-spectrum implementations;
- analytic Gaussian-beam, thin-lens, reversibility, conservation, and
  cross-model validation.

Exit gate:

- legacy quantum results remain within their recorded tolerances;
- analytic optics errors and trusted numerical regimes are reported;
- a CLI config can propagate a Gaussian beam through a lens or phase mask.

### M2 - Differentiable Solver v0.3

Deliverables:

- optional JAX dependencies and array-native operators;
- NumPy/JAX value parity;
- scan-based propagation suitable for JIT compilation;
- normalized mode-overlap and intensity objectives;
- centered directional finite-difference gradient tests;
- documented CPU and optional accelerator behavior.

Exit gate:

- required gradients meet declared relative-error tolerances;
- a compiled optimization step runs in the supported CPU development setup;
- GPU availability is not required for correctness or CI.

### M3 - Inverse-Design MVP v0.4

Deliverables:

- bounded phase-mask parameterization;
- smoothing, total-variation, and optional phase-quantization constraints;
- Adam-based optimization with deterministic seeds, restart support, early
  stopping, and best-design checkpoints;
- optimization history, resolved config, environment, fields, metrics, and
  report artifacts;
- a small deterministic CI optimization that improves its objective.

Exit gate:

- optimization results can be reproduced from committed inputs;
- the objective, constraints, and gradient evidence are visible in the report;
- interruptions can resume from a checkpoint without changing the experiment.

### M4 - Robust Flagship v1.0

Flagship question:

> Can a low-cost differentiable scalar propagator design a multi-plane
> phase-only Gaussian-to-HG10 converter that remains useful under held-out
> wavelength, alignment, phase-depth, and plane-spacing perturbations?

Required comparisons:

- no phase masks;
- one optimized phase mask;
- three masks optimized for nominal conditions;
- three masks optimized over a robustness ensemble.

Validation ladder:

1. analytic optics checks;
2. NumPy/JAX agreement;
3. grid, padding, aperture, and step convergence;
4. Fresnel-designed device evaluated with band-limited angular spectrum;
5. an optional reduced-scale Maxwell/FDTD check only if computationally
   practical.

Exit gate:

- optimization and held-out perturbation sets are distinct;
- nominal and worst-case performance are reported against every baseline;
- model-transfer loss and the trusted regime are explicit;
- CLI, dashboard, and README lead with the wave inverse-design identity while
  preserving `quantum-lab` compatibility.

## Initial Epic Issues

Create these as GitHub tracking issues and link child implementation issues from
their checklists:

1. `[Epic] Freeze and document the quantum v0.1 baseline`
2. `[Epic] Extract the functional wave-propagation core`
3. `[Epic] Add validated scalar-optics propagation`
4. `[Epic] Add JAX differentiation and gradient validation`
5. `[Epic] Implement the constrained inverse-design engine`
6. `[Epic] Build the robust multi-plane mode converter`
7. `[Epic] Validate transfer between propagation models`
8. `[Epic] Reorganize the CLI, dashboard, and documentation`

## Post-M0 Implementation Queue

After M0 is merged, create and execute narrowly scoped issues in this order:

1. Define the backend-neutral propagation interface.
2. Refactor the current solver to delegate to the functional kernel.
3. Add arbitrary complex input fields.
4. Add optical grids and physical-unit configuration.
5. Add Gaussian and Hermite-Gaussian sources.
6. Implement and validate Fresnel propagation.
7. Implement and validate band-limited angular-spectrum propagation.
8. Add optional JAX execution and NumPy/JAX parity tests.

## Tracking Rules

- GitHub issues are the source of truth for active status and ownership.
- Milestones group issues into releasable capabilities.
- One issue should normally map to one branch and one pull request.
- Architecture changes require an ADR or a superseding ADR.
- Numerical or scientific issues must state their validation evidence before
  implementation begins.
- Large outputs remain ignored. Commit configs, regeneration code, and small
  reference metrics; attach or upload full reports as workflow artifacts.
- Update this roadmap only when direction, milestone scope, status, or exit gates
  change.
