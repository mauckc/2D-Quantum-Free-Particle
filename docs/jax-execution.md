# JAX Execution and Precision

NumPy remains the correctness reference. JAX is an optional array-native path
for automatic differentiation and compilation; importing or running the legacy
quantum and NumPy optics workflows does not require it.

Install the CPU development path with:

```bash
uv sync --dev --extra jax
```

The project correctness tests run in 64-bit mode:

```bash
JAX_ENABLE_X64=1 uv run --extra jax pytest
```

On PowerShell, set `$env:JAX_ENABLE_X64 = "1"` before the command. JAX defaults
to 32-bit types, so callers that require reference comparisons must enable x64
before creating arrays or compiling functions. See the official JAX
[default-dtype guidance](https://docs.jax.dev/en/latest/default_dtypes.html) and
[installation matrix](https://docs.jax.dev/en/latest/installation.html).

`multiplane_scan_jax` uses `jax.lax.scan` over fixed-shape phase/propagation
pairs. It is suitable for `jax.jit` and later reverse-mode differentiation. The
CPU path is the supported correctness and CI configuration. An accelerator may
improve throughput, but GPU availability and GPU parity are not release gates.

The current JAX module implements Fresnel propagation. Band-limited angular
spectrum remains an independently formulated NumPy evaluation model so later
model-transfer checks do not reuse the optimizing implementation.
