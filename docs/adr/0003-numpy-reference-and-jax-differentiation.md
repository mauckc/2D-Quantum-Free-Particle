# ADR 0003: Keep NumPy as Reference and Add JAX for Differentiation

Status: Accepted

Date: 2026-07-13

## Context

The current solver and validation suite use NumPy, with optional FFT backends.
Differentiation requires every traced array operation—not only the FFT—to remain
inside the differentiable array system. Replacing the reference path outright
would make it harder to separate modeling errors from autodiff or compilation
errors.

## Decision

NumPy remains the correctness reference. JAX will be an optional differentiable
execution path using pure functions and array-native operators. Correctness and
CPU CI will not require a GPU. Accelerator execution is an optional performance
path.

## Consequences

- The existing FFT-only backend abstraction is insufficient for the new core;
  exponentials, reductions, casts, masks, and FFTs must all be backend-native.
- NumPy/JAX parity and directional finite-difference checks are required.
- 64-bit JAX execution is used for correctness comparisons where appropriate;
  lower precision may be offered for production optimization.
- CuPy and pyFFTW may continue serving forward simulations but do not
  automatically become differentiable backends.

## Alternatives Considered

- JAX-only implementation: rejected because it would remove the independent
  reference path.
- PyTorch as the first differentiable backend: viable, but not selected because
  JAX's NumPy-like functional transformations align more directly with the
  planned propagator and robustness ensembles.
- Hand-derived adjoints: deferred until profiling shows reverse-mode autodiff is
  inadequate for the target problem sizes.
