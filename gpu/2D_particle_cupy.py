"""GPU-accelerated 2D quantum free particle simulation using CuPy.

This script mirrors the Python implementation but executes all heavy
array operations on a CUDA-capable GPU. It should run on modern NVIDIA
GPUs such as the RTX 5090 as long as CuPy is installed with the
appropriate CUDA toolkit.
"""

import cupy as cp
import numpy as np
import math
import matplotlib.pyplot as plt


def compute_double_slit_potential(X, Y, L, barrier_width, slit_width,
                                  slit_distance, barrier_potential):
    """Return the double-slit potential evaluated on arrays X and Y."""
    barrier_center = L / 2
    within_barrier = (X > barrier_center - barrier_width / 2) & (
        X < barrier_center + barrier_width / 2
    )
    half_slit = slit_width / 2
    half_dist = slit_distance / 2
    in_first = (Y > barrier_center - half_slit - half_dist) & (
        Y < barrier_center - half_slit + half_dist
    )
    in_second = (Y > barrier_center + half_slit - half_dist) & (
        Y < barrier_center + half_slit + half_dist
    )
    in_slit = in_first | in_second
    return cp.where(within_barrier & (~in_slit), barrier_potential, 0.0)


def main():
    N = 128
    L = 10.0
    dt = 0.001
    tf = 5.0
    dx = L / N

    barrier_width = 0.5
    slit_width = 1.0
    slit_distance = 0.5
    barrier_potential = 200.0
    kx_magnitude = -5.0
    ky_magnitude = 0.0
    amplitude = 1.0
    sigma = 0.73

    x = cp.linspace(-L / 2, L / 2, N)
    y = cp.linspace(-L / 2, L / 2, N)
    X, Y = cp.meshgrid(x, y)

    V = compute_double_slit_potential(
        X, Y, L, barrier_width, slit_width, slit_distance, barrier_potential
    )

    kx = (kx_magnitude * math.pi) / L
    ky = (ky_magnitude * math.pi) / L
    gaussian = cp.exp(-((X * X) + (Y * Y)) / (4 * sigma * sigma))
    psi_real = amplitude * cp.cos(kx * X + ky * Y) * gaussian
    psi_imag = amplitude * cp.sin(kx * X + ky * Y) * gaussian

    k_vals = cp.fft.fftfreq(N, d=dx) * 2 * cp.pi
    kx_vals = k_vals.reshape(1, N)
    ky_vals = k_vals.reshape(N, 1)
    ksq = kx_vals ** 2 + ky_vals ** 2

    t = 0.0
    while t < tf:
        psi_real, psi_imag = (
            psi_real * cp.cos(V * dt) + psi_imag * cp.sin(V * dt),
            psi_imag * cp.cos(V * dt) - psi_real * cp.sin(V * dt),
        )

        psi_complex = psi_real + 1j * psi_imag
        chi_complex = cp.fft.fft2(psi_complex)

        phase = cp.exp(-1j * ksq * dt / 2)
        chi_complex *= phase

        psi_complex = cp.fft.ifft2(chi_complex)
        psi_real = cp.real(psi_complex)
        psi_imag = cp.imag(psi_complex)

        t += dt

    prob = psi_real ** 2 + psi_imag ** 2
    prob_cpu = cp.asnumpy(prob)
    plt.imshow(prob_cpu, cmap="hot", interpolation="nearest")
    plt.title("Final probability density")
    plt.show()


if __name__ == "__main__":
    main()
