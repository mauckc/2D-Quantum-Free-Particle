# v1.0 Release Evidence

The v1.0 release completes M0 through M4 without removing the legacy quantum
surface. Its release claim is deliberately narrow: a differentiable
monochromatic coherent scalar Fresnel model designed a three-plane phase-only
Gaussian-to-HG10 converter whose held-out worst overlap exceeds the nominally
optimized three-plane design, and the fixed device transfers to an independently
implemented band-limited angular-spectrum evaluator inside declared numerical
limits.

## Gate summary

- M0: quantum v0.1 references preserve free, barrier, measured-barrier, and
  double-slit behavior without treating their current diagnostics as new
  physical validation.
- M1: the pure propagation core, quantum adapter, scalar sources/elements,
  analytic Gaussian/lens checks, conservation/reversibility, convergence, and
  Fresnel/BLAS reference regime pass.
- M2: JAX x64 parity, fixed-shape JIT scan, normalized objectives, centered
  directional gradient checks, and a compiled CPU update pass.
- M3: bounded/smoothed/TV/quantized constraints, deterministic Adam improvement,
  early stopping, best checkpoints, fingerprinted restart equivalence, and
  complete artifacts pass.
- M4: all four baselines, disjoint robust/held-out sets, worst-case improvement,
  BLAS transfer, retained spectrum, grid/padding/aperture/substep convergence,
  boundary power, CLI, dashboard, documentation, and CI pass.

The small committed references live in `benchmarks/reference/`; regeneration
code is adjacent. Generated fields, histories, checkpoints, and reports remain
ignored under `runs/` and `reports/`.

## Compatibility and limitations

`wave-lab` is the primary v1.0 CLI. `quantum-lab` remains available with the
same commands, committed configs, artifact schema, and dashboard workflows.
The Python import namespace and distribution name remain unchanged to avoid a
compatibility-breaking rename in this release.

No vector-Maxwell, polarization, broadband, nonlinear, partially coherent,
high-index integrated-photonics, or manufacturing-readiness claim is made.
Maxwell/FDTD remains optional and was not computationally proportionate for the
millimetre-scale CPU-CI flagship; the accepted release gate is independent
scalar-model transfer under ADR 0004.
