"""Regenerate the small robust-flagship M4 comparison reference."""

from __future__ import annotations

import json
from pathlib import Path

from quantum_dynamics_lab.flagship import run_flagship
from quantum_dynamics_lab.flagship_config import load_flagship_config


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPOSITORY_ROOT / "examples" / "robust_flagship.toml"
REFERENCE_PATH = Path(__file__).with_name("flagship_m4.json")


def generate_metrics() -> dict:
    return run_flagship(load_flagship_config(CONFIG_PATH)).metrics


def main() -> None:
    REFERENCE_PATH.write_text(
        json.dumps(generate_metrics(), indent=2) + "\n",
        encoding="utf-8",
    )
    print(REFERENCE_PATH)


if __name__ == "__main__":
    main()
