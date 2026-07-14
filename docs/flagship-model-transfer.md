# Flagship Model Transfer and Trusted Regime

The three phase masks are optimized only with differentiable Fresnel
propagation. The committed validation re-evaluates every required baseline under
nominal and every held-out scenario with the independently formulated NumPy
band-limited angular-spectrum (BLAS) model. The code does not reuse the JAX
Fresnel transfer function for that evaluation.

The report records per-scenario overlap and input-referenced efficiency from
both models, absolute transfer loss, and the minimum spectral power retained by
BLAS. It also fixes the robust design and runs:

- 24, 32, 48, and 64 sample grid studies at fixed window size;
- 1.0×, 1.5×, and 2.0× padding at fixed sample spacing;
- 1.2, 1.5, 1.8, and 1.95 mm aperture radii;
- one, two, and four propagation substeps per inter-mask distance;
- a padded-window boundary-power diagnostic.

The thresholds are declared in `flagship_validation.DECLARED_LIMITS` and stored
with every generated report. A result is trusted only if all transfer,
retained-spectrum, convergence, and boundary checks pass. Regenerate the small
reference with:

```bash
JAX_ENABLE_X64=1 uv run --extra jax python -m \
  benchmarks.reference.generate_flagship_validation_m4
```

The trusted regime remains monochromatic coherent scalar free-space optics at
the committed aperture, sampling, window, and angular content. Reduced-scale
Maxwell/FDTD is not performed: the millimetre-scale device is outside practical
CPU CI scope, and ADR 0004 makes independent scalar-model transfer—not optional
full-wave simulation—the v1.0 release gate.
