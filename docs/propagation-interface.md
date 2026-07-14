# Backend-Neutral Propagation Interface

Issue #27 introduces a pure symmetric split-step interface in
`quantum_dynamics_lab.propagation`. It does not yet change the maintained
`SplitStepSolver`; that delegation belongs to issue #28 so parity can be
measured against the pre-refactor implementation.

## Responsibilities

The caller owns:

- the physical model, units, grid, and evolution step;
- construction of position-space and spectral generator arrays;
- boundary masks, normalization, diagnostics, and experiment branching;
- configuration parsing, plotting, artifact writing, and dashboards.

The propagation module owns:

- backend-native exponentiation of the generator arrays;
- the half-position, full-spectral, half-position composition;
- validation that a field and its operators are matching 2D arrays.

The explicit `PropagationBackend` owns array conversion, dtype/device policy,
elementwise operations, and FFT execution. `NumPyPropagationBackend` is the
complex128 correctness reference. Future JAX support must implement the same
operations without converting traced arrays through NumPy.

## Generator Convention

For an evolution equation written as

```text
d(field) / ds = (G_position + G_spectral) field,
```

`prepare_split_step_operators` constructs

```text
exp(0.5 * ds * G_position)
exp(      ds * G_spectral)
exp(0.5 * ds * G_position)
```

The generators and `step_size` may use any consistent physical units. A
negative step size prepares the inverse phases for reversibility checks.

## NumPy Reference Example

For the existing dimensionless quantum model, an adapter can prepare generators
without passing configuration objects into the kernel:

```python
import numpy as np

from quantum_dynamics_lab.propagation import (
    NumPyPropagationBackend,
    prepare_split_step_operators,
    propagate_split_step,
)

backend = NumPyPropagationBackend()
position_generator = -1j * potential / hbar
spectral_generator = -1j * hbar * k_squared / (2.0 * mass)
operators = prepare_split_step_operators(
    position_generator,
    spectral_generator,
    dt,
    backend=backend,
)
next_field = propagate_split_step(field, operators, backend=backend)
```

Boundary damping and absorbed-probability diagnostics are intentionally applied
by the quantum adapter after `next_field` is returned.

## Validation Contract

`tests/test_propagation.py` records the interface-level evidence required before
the maintained solver delegates to this module:

- analytic phase evolution for a periodic plane wave;
- norm conservation and forward/reverse reconstruction for a real potential;
- decreasing free-Gaussian variance error under grid refinement;
- complex-array parity with the untouched solver for free and static potentials;
- matching 2D shape enforcement and non-mutation of input arrays.

The committed quantum v0.1 metrics remain the end-to-end compatibility gate.
