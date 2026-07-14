# Robust Gaussian-to-HG10 Flagship

The committed `examples/robust_flagship.toml` asks whether a low-cost
differentiable scalar propagator can design a phase-only Gaussian-to-HG10 mode
converter that remains useful under held-out wavelength, alignment, phase-depth,
and plane-spacing perturbations.

The required comparison uses common source, target, aperture, total length,
sampling, constraints, metrics, and seed:

1. no phase masks;
2. one optimized phase mask;
3. three phase masks optimized nominally;
4. three phase masks optimized over the declared robustness ensemble.

The optimization ensemble includes nominal conditions plus distinct variations
of all four perturbation classes. The held-out set uses different values and, for
alignment, diagonal directions not used during optimization. Scenario signatures
are asserted disjoint in tests and recorded with every generated report.

Run the comparison with:

```bash
JAX_ENABLE_X64=1 uv run --extra jax wave-lab flagship \
  examples/robust_flagship.toml --out runs/robust-flagship
```

The initial report contains nominal, held-out mean, and held-out worst overlap,
intensity similarity, and input-referenced mode efficiency for all four
baselines. It is still a Fresnel-only result. The M4 conclusion is not complete
until independent BLAS model transfer and grid, padding, aperture, and
propagation-step convergence are recorded.
