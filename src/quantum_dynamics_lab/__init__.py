"""Differentiable 2D wave inverse-design lab with preserved quantum demos."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("quantum-dynamics-lab")
except PackageNotFoundError:  # pragma: no cover - editable tree fallback
    __version__ = "1.0.0"
