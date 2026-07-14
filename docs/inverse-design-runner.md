# M3 Inverse-Design Runner

The M3 runner optimizes one phase mask per declared propagation distance using
JAX x64 Fresnel propagation and deterministic Adam. Unconstrained variables are
smoothed, bounded, and optionally straight-through quantized during training;
saved and reported designs always use hard quantization. The loss can combine
mode-overlap, intensity-similarity, and wrapped-TV terms.

Run the committed deterministic smoke experiment with:

```bash
JAX_ENABLE_X64=1 uv run --extra jax wave-lab optimize \
  examples/inverse_design_smoke.toml --out runs/inverse-design-smoke
```

Resume an interrupted run without changing its config:

```bash
JAX_ENABLE_X64=1 uv run --extra jax wave-lab optimize \
  examples/inverse_design_smoke.toml --out runs/inverse-design-smoke \
  --resume runs/inverse-design-smoke/checkpoint.npz
```

The checkpoint contains parameters, both Adam moments, iteration, best design,
early-stop state, and full history. A SHA-256 fingerprint of the complete
resolved config prevents accidental restart under changed experiment semantics.

Each run writes `history.json`, `resolved_config.json`, `environment.json`,
`fields.npz`, `metrics.json`, `best_design.npz`, `checkpoint.npz`, and
`report.md`. Provenance includes git commit/dirty state, dependency versions,
backend/device, x64 status, and seed. The report exposes objectives, constraints,
and the committed M2 gradient evidence.

The deterministic CI smoke establishes objective improvement and restart
equivalence for the discrete optimizer. It does not establish robustness or
physical transfer; those are M4 held-out and independent-model requirements.
