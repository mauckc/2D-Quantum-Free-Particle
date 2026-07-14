# Phase Parameterization and Constraints

Unconstrained real design variables are mapped to phase with a shifted/scaled
`tanh`, keeping every value inside explicit bounds. Optional periodic five-point
smoothing is applied to the design variables before the bounded map, so the
result remains bounded while high-frequency roughness is reduced.

Wrapped isotropic total variation uses nearest-neighbor phase differences mapped
through `atan2(sin(delta), cos(delta))`. It therefore measures physical phase
jumps modulo 2π and is invariant when integer phase cycles are added to samples.
It is a regularization term, not a hard fabrication guarantee.

Optional hard quantization maps phase to equally spaced levels in the half-open
interval `[phase_min, phase_max)`. Evaluation uses the hard values. Optimization
may use a straight-through estimator: its forward value is exactly quantized,
while its local backward derivative is the identity. Reports must disclose when
this biased estimator is used and must evaluate the saved design with hard
quantization.

NumPy supplies the reference functions in `constraints.py`; JAX supplies the
JIT-compatible differentiable path in `jax_constraints.py`. Focused tests cover
array parity, exact levels, bounds, roughness reduction, wrapped-TV invariance,
straight-through gradients, and a centered directional derivative through the
composed constraint path.
