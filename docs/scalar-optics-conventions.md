# Scalar-Optics Conventions

The optics layer models a monochromatic, coherent scalar field sampled on a
transverse periodic grid. It is a free-space/paraxial laboratory model, not a
vector-Maxwell or high-index integrated-photonics solver.

Lengths, wavelength, waist, and focal length may use any one consistent unit;
the committed examples use SI units. Arrays are indexed as `(y, x)`. Spatial
frequencies are cycles per unit length in NumPy FFT order. Integrated power is

```text
sum(abs(E)**2) * dx * dy.
```

Gaussian and Hermite-Gaussian constructors return numerically normalized
complex128 fields. Hermite-Gaussian orders use the conventional `(x, y)` order,
so HG10 is `(1, 0)`. The Gaussian `waist` is the field-amplitude 1/e radius.
Arbitrary fields are copied without implicit normalization so measured or
externally generated amplitudes retain their meaning.

For forward propagation with carrier convention `exp(+i k z)`, a thin lens has
transmission

```text
exp(-i * k * ((x-x0)^2 + (y-y0)^2) / (2*f)).
```

Apertures are passive binary transmissions. Phase-only masks must lie inside
their declared phase bounds and have unit magnitude. Thin-element application
rejects gain, does not mutate its inputs, and remains separate from propagation,
configuration, plotting, artifact writing, experiments, and dashboards.

## Forward CLI

Run the committed Gaussian-through-lens example with:

```bash
uv run wave-lab propagate examples/gaussian_lens.toml --out runs/gaussian-lens
```

The optics TOML schema is separate from the legacy quantum schema. Ordered
`[[steps]]` may apply lenses, circular or rectangular apertures, bounded
sinusoidal phase masks, and Fresnel or BLAS propagation. The output directory
contains input/final fields, resolved config, provenance, and concise metrics.
The existing `quantum-lab` commands and artifacts remain available unchanged.
