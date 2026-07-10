"""2D Quantum Dynamics Lab."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("quantum-dynamics-lab")
except PackageNotFoundError:  # pragma: no cover - editable tree fallback
    __version__ = "0.1.0"
