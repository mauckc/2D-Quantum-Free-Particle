from __future__ import annotations

import tomllib
from pathlib import Path

from quantum_dynamics_lab import __version__


ROOT = Path(__file__).resolve().parents[1]


def test_v1_package_and_public_entrypoints_are_declared() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]

    assert project["version"] == "1.0.0"
    assert __version__ == "1.0.0"
    assert project["scripts"] == {
        "quantum-lab": "quantum_dynamics_lab.cli:main",
        "wave-lab": "quantum_dynamics_lab.wave_cli:main",
    }


def test_m4_release_identity_and_evidence_are_durable() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    roadmap = (ROOT / "ROADMAP.md").read_text(encoding="utf-8")

    assert readme.startswith("# Differentiable 2D Wave Inverse-Design Lab")
    assert "`quantum-lab`" in readme
    assert "| M4 - Robust flagship v1.0 | Complete |" in roadmap
    assert (ROOT / "benchmarks/reference/flagship_m4.json").is_file()
    assert (ROOT / "benchmarks/reference/flagship_validation_m4.json").is_file()
    assert (ROOT / "docs/v1-release.md").is_file()
