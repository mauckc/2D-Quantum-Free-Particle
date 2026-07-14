# Optical Objectives and Gradient Evidence

The objective layer exposes three dimensionless metrics:

- **Mode overlap** is the squared magnitude of the complex inner product divided
  by input and target powers. It is insensitive to global phase, equals one for
  identical modes, and equals zero for orthogonal modes.
- **Intensity similarity** is cosine similarity between sampled intensities. It
  measures shape agreement without requiring phase agreement.
- **Mode power efficiency** is power coupled into the declared target mode
  divided by the incident input power. Unlike normalized overlap, it penalizes
  loss before the output plane.

For zero-power inputs or targets the metric is undefined and returns NaN. NumPy
validates matching finite 2D arrays and positive sample area before evaluation.
JAX uses array-native operations and produces NaN inside compiled code for the
same undefined cases.

## M2 Gradient Check

`benchmarks/reference/jax_m2.json` records a deterministic three-plane,
32-by-32 Gaussian-to-HG10 case with seed `20260714`. For both mode-overlap and
intensity losses, it compares the reverse-mode directional derivative with

```text
(loss(theta + epsilon*d) - loss(theta - epsilon*d)) / (2*epsilon)
```

using a normalized deterministic direction and `epsilon = 1e-5` in JAX x64.
The declared limits are `1e-9` absolute and `1e-6` relative error. The evidence
also records NumPy/JAX objective parity and executes one JIT-compiled finite
parameter update on CPU.

Regenerate the evidence after installing the JAX extra:

```bash
JAX_ENABLE_X64=1 uv run --extra jax python -m benchmarks.reference.generate_jax_m2
```

This validates derivatives of the implemented discrete Fresnel model. It does
not validate the scalar physical approximation or establish optimizer
robustness; those require the separate model-transfer and held-out ladder.
