#include <cuda_runtime.h>
#include <cufft.h>
#include <stdio.h>
#include <math.h>

// Simple GPU-based 2D quantum free particle solver
// Uses a split-step scheme with cuFFT. Parameters are
// chosen for demonstration and can be adjusted.

#define N 128
#define L 10.0
#define DT 0.001
#define TF 5.0
#define SIGMA 0.73
#define AMPLITUDE 1.0

#define BARRIER_WIDTH 0.5
#define SLIT_WIDTH 1.0
#define SLIT_DISTANCE 0.5
#define BARRIER_POTENTIAL 200.0

__global__ void setup_grid(double *x, double *y)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        double dx = L / N;
        x[idx] = -L / 2.0 + idx * dx;
        y[idx] = -L / 2.0 + idx * dx;
    }
}

__global__ void init_potential(const double *x, const double *y, double *V)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < N && j < N)
    {
        double barrier_center = L / 2.0;
        double X = x[i];
        double Y = y[j];

        bool within_barrier = (X > barrier_center - BARRIER_WIDTH / 2.0) &&
                              (X < barrier_center + BARRIER_WIDTH / 2.0);
        double half_slit = SLIT_WIDTH / 2.0;
        double half_dist = SLIT_DISTANCE / 2.0;
        bool in_first = (Y > barrier_center - half_slit - half_dist) &&
                        (Y < barrier_center - half_slit + half_dist);
        bool in_second = (Y > barrier_center + half_slit - half_dist) &&
                         (Y < barrier_center + half_slit + half_dist);
        bool in_slit = in_first || in_second;
        V[j * N + i] = (within_barrier && !in_slit) ? BARRIER_POTENTIAL : 0.0;
    }
}

__global__ void init_wavefunction(const double *x, const double *y, cufftDoubleComplex *psi)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < N && j < N)
    {
        double kx = (-5.0 * M_PI) / L;
        double ky = 0.0;
        double X = x[i];
        double Y = y[j];
        double gaussian = exp(-((X * X) + (Y * Y)) / (4.0 * SIGMA * SIGMA));
        double phase = kx * X + ky * Y;
        psi[j * N + i].x = AMPLITUDE * cos(phase) * gaussian;
        psi[j * N + i].y = AMPLITUDE * sin(phase) * gaussian;
    }
}

__global__ void update_position(cufftDoubleComplex *psi, const double *V)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < N && j < N)
    {
        double phase = V[j * N + i] * DT;
        double c = cos(phase);
        double s = sin(phase);
        double re = psi[j * N + i].x;
        double im = psi[j * N + i].y;
        psi[j * N + i].x = re * c + im * s;
        psi[j * N + i].y = im * c - re * s;
    }
}

__global__ void update_momentum(cufftDoubleComplex *chi, const double *kx, const double *ky)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < N && j < N)
    {
        double ksq = kx[i] * kx[i] + ky[j] * ky[j];
        double phase = (ksq * DT) / 2.0;
        double c = cos(phase);
        double s = sin(phase);
        double re = chi[j * N + i].x;
        double im = chi[j * N + i].y;
        chi[j * N + i].x = im * s + re * c;
        chi[j * N + i].y = im * c - re * s;
    }
}

__global__ void normalize(cufftDoubleComplex *psi)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N * N)
    {
        psi[idx].x /= (N * N);
        psi[idx].y /= (N * N);
    }
}

int main()
{
    double *x, *y, *V, *kx, *ky;
    cufftDoubleComplex *psi, *chi;

    cudaMallocManaged(&x, N * sizeof(double));
    cudaMallocManaged(&y, N * sizeof(double));
    cudaMallocManaged(&V, N * N * sizeof(double));
    cudaMallocManaged(&kx, N * sizeof(double));
    cudaMallocManaged(&ky, N * sizeof(double));
    cudaMallocManaged(&psi, N * N * sizeof(cufftDoubleComplex));
    cudaMallocManaged(&chi, N * N * sizeof(cufftDoubleComplex));

    setup_grid<<<(N + 255) / 256, 256>>>(x, y);
    cudaDeviceSynchronize();

    for (int i = 0; i < N; ++i)
    {
        kx[i] = 2.0 * M_PI * (((i + N / 2) % N) - N / 2) / L;
        ky[i] = 2.0 * M_PI * (((i + N / 2) % N) - N / 2) / L;
    }

    dim3 blocks((N + 15) / 16, (N + 15) / 16);
    dim3 threads(16, 16);
    init_potential<<<blocks, threads>>>(x, y, V);
    init_wavefunction<<<blocks, threads>>>(x, y, psi);
    cudaDeviceSynchronize();

    cufftHandle plan;
    cufftPlan2d(&plan, N, N, CUFFT_Z2Z);

    double t = 0.0;
    while (t < TF)
    {
        update_position<<<blocks, threads>>>(psi, V);
        cudaDeviceSynchronize();

        cufftExecZ2Z(plan, psi, chi, CUFFT_FORWARD);
        cudaDeviceSynchronize();

        update_momentum<<<blocks, threads>>>(chi, kx, ky);
        cudaDeviceSynchronize();

        cufftExecZ2Z(plan, chi, psi, CUFFT_INVERSE);
        cudaDeviceSynchronize();

        normalize<<<(N * N + 255) / 256, 256>>>(psi);
        cudaDeviceSynchronize();

        t += DT;
    }

    int center = (N / 2) * N + (N / 2);
    double prob = psi[center].x * psi[center].x + psi[center].y * psi[center].y;
    printf("Probability at center: %f\n", prob);

    cufftDestroy(plan);
    cudaFree(chi);
    cudaFree(psi);
    cudaFree(V);
    cudaFree(x);
    cudaFree(y);
    cudaFree(kx);
    cudaFree(ky);
    return 0;
}

