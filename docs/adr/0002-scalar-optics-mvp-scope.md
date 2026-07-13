# ADR 0002: Use Scalar Free-Space Optics for the MVP

Status: Accepted

Date: 2026-07-13

## Context

The maintained solver evolves a complex scalar field with FFT-based kinetic or
diffraction operators and pointwise phase operators. Scalar free-space optics
and thin phase elements reuse that structure directly. High-index integrated
photonics would require vector fields, polarization, open ports, dispersive
materials, and higher-fidelity Maxwell solvers.

## Decision

The first photonics application will be monochromatic, coherent, scalar
free-space propagation through apertures, lenses, phase masks, and optional
smooth index variations. The first flagship design will be a multi-plane
phase-only Gaussian-to-HG10 mode converter.

## Consequences

- Fresnel/paraxial and band-limited angular-spectrum models are in scope.
- Polarization, vector Maxwell fields, nonlinear optics, and manufacturing-ready
  layout are out of scope for the MVP.
- Reports must state the scalar approximation and its validated regime.
- Integrated-photonics examples may be added later only with suitable validation.

## Alternatives Considered

- Semiconductor quantum transport: rejected for the first pivot because it
  requires contacts, electrostatics, and NEGF or comparable transport models.
- High-index silicon-photonic inverse design: deferred because it is not a
  trustworthy application of the present scalar propagator.
- Continuing with Quantum Zeno research: retained only as legacy demonstration
  because credible extension would require open-system measurement models.
