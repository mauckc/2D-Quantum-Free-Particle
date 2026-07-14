# ADR 0001: Preserve Quantum Compatibility During the Pivot

Status: Accepted

Date: 2026-07-13

## Context

The repository has a working Python package, CLI, dashboard, experiment configs,
reports, and regression suite for 2D quantum demonstrations. Replacing them in a
single rewrite would discard useful numerical coverage and make regressions hard
to identify.

## Decision

The wave inverse-design functionality will be added incrementally. Existing
quantum experiments will become adapters and legacy demonstrations over the new
generic propagation core. The `quantum-lab` command and committed configs remain
available through the initial optics and inverse-design releases.

## Consequences

- Core refactors require parity tests against recorded quantum results.
- Public documentation must distinguish legacy demonstrations from the new
  research direction.
- Some compatibility code will remain temporarily after the new `wave-lab`
  surface is introduced.
- Package renaming is deferred until the new core and flagship experiment work.

## Alternatives Considered

- A clean-slate optics repository: rejected because it would strand tested code
  and split maintenance.
- Immediate package and CLI renaming: rejected because it would combine branding
  churn with numerical refactoring.
