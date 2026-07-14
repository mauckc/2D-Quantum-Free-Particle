# Reference Metrics

This directory holds small, committed numerical reference data used to detect
regressions and document validated ranges.

Reference files should contain metrics and tolerances, not large wave fields,
videos, or generated reports. Each file must identify:

- the generating config;
- the metric definition;
- the expected value or acceptable range;
- the numerical tolerance and its justification;
- the implementation or milestone that established the reference.

Large run outputs belong under ignored `runs/` or `reports/` directories and may
be published as CI artifacts or release assets.

## Quantum v0.1 baseline

`quantum_v0_1.json` freezes scalar metrics from the committed free-packet,
barrier/Zeno, and double-slit example configurations. The reference generator
forces the NumPy complex128 backend and includes both unmeasured and measured
variants of the barrier configuration.

Regenerate the file from the repository root with:

```bash
uv run python benchmarks/reference/generate_quantum_v0_1.py
```

The recorded tolerances allow only floating-point/FFT roundoff across supported
platforms. They are regression tolerances, not uncertainty estimates and not
evidence that the current Zeno or double-slit interpretation is physically
validated.

## Scalar optics M1

`optics_m1.json` records the analytic Gaussian, thin-lens focus, conservation,
reversibility, grid/window convergence, retained-band, and Fresnel-to-angular-
spectrum metrics that establish the M1 NumPy reference regime.

Regenerate the file from the repository root with:

```bash
uv run python -m benchmarks.reference.generate_optics_m1
```

The limits are numerical acceptance thresholds for coherent monochromatic
scalar free-space optics. Their trusted regime is documented in
`docs/optics-model-validation.md`; they do not establish vector-Maxwell or
high-index integrated-photonics accuracy.

## Differentiable solver M2

`jax_m2.json` records NumPy/JAX objective parity, centered directional finite-
difference checks for mode and intensity losses, and a compiled CPU update. It
uses deterministic seed `20260714` and JAX x64. Regenerate it with:

```bash
JAX_ENABLE_X64=1 uv run --extra jax python -m benchmarks.reference.generate_jax_m2
```

The gradient tolerances validate the implemented discrete Fresnel computation,
not the scalar physical model or an optimization result.
