"""NumPy reference objectives for scalar optical fields."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike


__all__ = [
    "intensity_similarity",
    "mode_overlap",
    "mode_power_efficiency",
]


def mode_overlap(
    field: ArrayLike,
    target: ArrayLike,
    sample_area: float = 1.0,
) -> float:
    """Return normalized phase-insensitive complex-mode overlap in [0, 1]."""

    values, target_values, area = _validated_fields(field, target, sample_area)
    inner = np.vdot(target_values, values) * area
    field_power = float(np.sum(np.abs(values) ** 2) * area)
    target_power = float(np.sum(np.abs(target_values) ** 2) * area)
    denominator = field_power * target_power
    if denominator <= 0.0:
        return float("nan")
    return float(np.abs(inner) ** 2 / denominator)


def intensity_similarity(
    field: ArrayLike,
    target: ArrayLike,
    sample_area: float = 1.0,
) -> float:
    """Return cosine similarity of field intensities in [0, 1]."""

    values, target_values, area = _validated_fields(field, target, sample_area)
    intensity = np.abs(values) ** 2
    target_intensity = np.abs(target_values) ** 2
    numerator = float(np.sum(intensity * target_intensity) * area)
    denominator = math.sqrt(
        float(np.sum(intensity**2) * area)
        * float(np.sum(target_intensity**2) * area)
    )
    if denominator <= 0.0:
        return float("nan")
    return numerator / denominator


def mode_power_efficiency(
    field: ArrayLike,
    target: ArrayLike,
    input_power: float,
    sample_area: float = 1.0,
) -> float:
    """Return target-coupled output power divided by declared input power."""

    values, target_values, area = _validated_fields(field, target, sample_area)
    if not math.isfinite(input_power) or input_power <= 0.0:
        return float("nan")
    target_power = float(np.sum(np.abs(target_values) ** 2) * area)
    if target_power <= 0.0:
        return float("nan")
    inner = np.vdot(target_values, values) * area
    coupled_power = float(np.abs(inner) ** 2 / target_power)
    return coupled_power / input_power


def _validated_fields(
    field: ArrayLike,
    target: ArrayLike,
    sample_area: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    values = np.asarray(field, dtype=np.complex128)
    target_values = np.asarray(target, dtype=np.complex128)
    if values.ndim != 2 or values.shape != target_values.shape:
        raise ValueError("field and target must be matching two-dimensional arrays")
    area = float(sample_area)
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError("sample_area must be finite and positive")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(target_values)):
        raise ValueError("field and target must contain only finite values")
    return values, target_values, area
