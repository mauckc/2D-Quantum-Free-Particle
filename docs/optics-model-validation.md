# M1 Optics Model Validation

The NumPy reference layer contains two separately formulated free-space models:

- Fresnel propagation uses the paraxial transfer phase
  `exp(i*k*z - i*pi*lambda*z*(fx^2+fy^2))`.
- Band-limited angular spectrum (BLAS) uses exact scalar longitudinal dispersion
  `exp(i*k*z*sqrt(1-lambda^2*(fx^2+fy^2)))`, rejects evanescent samples, and
  applies distance/window-dependent alias-safe rectangular limits.

Both use periodic FFT boundaries. A field must be padded so that negligible
power reaches the transverse window boundary. Negative distance reverses phase;
BLAS reversibility additionally requires the input spectrum to lie inside the
retained band.

## Committed Evidence

`benchmarks/reference/optics_m1.json` records deterministic analytic Gaussian,
thin-lens focus, conservation, reversibility, grid/window convergence, retained
spectral power, and cross-model metrics. Regenerate it with:

```bash
uv run python -m benchmarks.reference.generate_optics_m1
```

The evidence uses a 1064 nm wavelength, sub-millimetre Gaussian waists,
12--18 cm propagation, and windows at least 6 mm across for the final
comparisons. Doubling both grid count and window width at fixed sampling reduces
the analytic Gaussian boundary error by more than four orders of magnitude.

## Trusted Regime

Treat a result as numerically trusted only when all of the following hold:

1. sampled input and output have negligible power at the window boundary;
2. BLAS retains effectively all significant spectral power (the M1 reference
   requires at least `1 - 1e-13`);
3. halving sample spacing and increasing padding do not materially change the
   reported metric;
4. Fresnel and BLAS agree at the tolerance appropriate to the claim (the M1
   paraxial reference requires relative complex-field and intensity errors below
   `5e-6`);
5. the physical problem is monochromatic, coherent, scalar, and free-space.

Fresnel error grows with angular content because it replaces exact longitudinal
dispersion with its paraxial expansion. BLAS is the independent scalar transfer
model for later model-transfer evaluation; agreement between the two does not
establish vector-Maxwell accuracy, polarization behavior, or validity for
high-index integrated photonics.
