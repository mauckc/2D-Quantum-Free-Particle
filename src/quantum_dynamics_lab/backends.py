from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import NDArray


ComplexArray = NDArray[np.complex128]


class FFTBackend(Protocol):
    name: str

    def fft2(self, values: ComplexArray) -> ComplexArray:
        ...

    def ifft2(self, values: ComplexArray) -> ComplexArray:
        ...


@dataclass(frozen=True)
class NumPyBackend:
    name: str = "numpy"

    def fft2(self, values: ComplexArray) -> ComplexArray:
        return np.fft.fft2(values)

    def ifft2(self, values: ComplexArray) -> ComplexArray:
        return np.fft.ifft2(values)


@dataclass(frozen=True)
class PyFFTWBackend:
    name: str = "pyfftw"

    def __post_init__(self) -> None:
        try:
            import pyfftw.interfaces.numpy_fft as numpy_fft  # noqa: F401
        except ImportError as exc:
            raise RuntimeError("pyFFTW backend requested but pyfftw is not installed") from exc

    def fft2(self, values: ComplexArray) -> ComplexArray:
        import pyfftw.interfaces.numpy_fft as fftw_fft

        return fftw_fft.fft2(values)

    def ifft2(self, values: ComplexArray) -> ComplexArray:
        import pyfftw.interfaces.numpy_fft as fftw_fft

        return fftw_fft.ifft2(values)


@dataclass(frozen=True)
class CuPyBackend:
    name: str = "cupy"

    def __post_init__(self) -> None:
        try:
            import cupy as cp  # noqa: F401
        except ImportError as exc:
            raise RuntimeError("CuPy backend requested but cupy is not installed") from exc

    def fft2(self, values: ComplexArray) -> ComplexArray:
        import cupy as cp

        return cp.asnumpy(cp.fft.fft2(cp.asarray(values)))

    def ifft2(self, values: ComplexArray) -> ComplexArray:
        import cupy as cp

        return cp.asnumpy(cp.fft.ifft2(cp.asarray(values)))


def create_backend(name: str) -> FFTBackend:
    normalized = name.lower()
    if normalized == "numpy":
        return NumPyBackend()
    if normalized == "pyfftw":
        return PyFFTWBackend()
    if normalized == "cupy":
        return CuPyBackend()
    if normalized == "auto":
        try:
            return PyFFTWBackend()
        except RuntimeError:
            return NumPyBackend()
    raise ValueError(f"unsupported FFT backend: {name}")
