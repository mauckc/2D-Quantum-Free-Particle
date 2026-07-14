# ADR 0004: Validate with a Propagation-Model Ladder

Status: Accepted

Date: 2026-07-13

## Context

An internally consistent optimizer can exploit discretization and model errors.
Successful optimization under one propagation function is not evidence that a
physical device or even a second numerical model will reproduce the result.
Full-wave Maxwell validation of millimetre-scale free-space systems may also be
computationally impractical.

## Decision

Validation will progress through an explicit ladder:

1. analytic Gaussian, lens, conservation, and reversibility checks;
2. NumPy/JAX value agreement;
3. grid, padding, aperture, and propagation-step convergence;
4. Fresnel-designed devices evaluated with an independently implemented
   band-limited angular-spectrum model;
5. optional reduced-scale Maxwell/FDTD validation when practical.

Held-out perturbations are required for robustness claims.

## Consequences

- Model-transfer loss becomes a primary reported metric.
- The trusted numerical regime must be reported alongside successful results.
- Maxwell validation is not a release blocker for the first scalar-optics MVP,
  but scalar-model transfer is.
- Research issues must declare their optimization and held-out evaluation sets
  before final results are interpreted.

## Alternatives Considered

- Validate only against the optimizing model: rejected because it cannot reveal
  model exploitation.
- Require full Maxwell validation for every milestone: rejected as disproportionate
  and potentially infeasible for the initial free-space scale.
