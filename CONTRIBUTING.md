# Contributing

## Before Starting

Read `AGENTS.md`, `ROADMAP.md`, and any relevant decision records under
`docs/adr/`. Find or create a GitHub issue whose outcome and acceptance criteria
match the proposed change.

An issue is ready when it has:

- one observable outcome;
- explicit scope and non-goals;
- measurable acceptance criteria;
- the required validation or scientific evidence;
- known dependencies;
- expected generated artifacts, when applicable.

## Branch and Pull-Request Workflow

Use one primary issue per branch and pull request. Branches created by Codex use
the form:

```text
codex/<issue-number>-<short-description>
```

Open a draft pull request early for work that spans more than one development
session. Use `Closes #<issue>` in the PR body when merging should close the
issue.

## Definition of Done

A change is done when:

- the issue acceptance criteria are satisfied;
- focused tests and the full regression suite pass;
- relevant configs and user-facing documentation are updated;
- required analytic, gradient, backend, convergence, or held-out evidence is
  recorded;
- generated artifacts include enough provenance to reproduce the result;
- compatibility changes are explicit;
- the PR states what its scientific evidence does and does not establish.

## Validation Commands

Minimum regression check:

```bash
uv sync --dev
uv run pytest
```

Dashboard checks:

```bash
uv sync --dev --extra ui
uv run --extra ui streamlit run apps/streamlit_dashboard.py
```

Additional commands belong in the issue acceptance criteria as optics,
differentiation, optimization, and benchmark workflows are added.

## Scientific and Numerical Changes

Match evidence to the claim:

- A refactor needs regression parity, not a new physical interpretation.
- A propagation model needs analytic checks and convergence evidence.
- A differentiable operator needs a directional finite-difference check.
- An optimizer needs deterministic improvement and complete history.
- A robustness claim needs held-out perturbations.
- A transfer claim needs an independently implemented evaluation model.

Do not weaken a tolerance solely to make a failing result pass. Explain and
justify any tolerance change in the issue and PR.

## Generated Artifacts

`runs/` and `reports/` are intentionally ignored. Do not commit large NPZ files,
videos, or transient optimizer checkpoints.

Commit instead:

- the complete input configuration;
- deterministic seeds;
- regeneration code;
- small reference metrics or expected ranges;
- concise documentation of the result.

Optimization and research runs should record the Git commit, dirty-state flag,
resolved config, dependency versions, backend/device, random seed, optimization
history, and final plus held-out metrics.
