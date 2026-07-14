# Differentiable 2D Wave Inverse-Design Lab

This repository is a reproducible 2D scalar-wave propagation and inverse-design
laboratory. NumPy is the correctness reference; optional JAX supplies automatic
differentiation and JIT-compiled multi-plane optimization. The v1.0 flagship is
a robust three-plane phase-only Gaussian-to-HG10 converter evaluated under
held-out perturbations and independently cross-checked with band-limited angular
spectrum propagation.

The original `quantum-lab` CLI, configs, dashboard experiments, and frozen
quantum v0.1 regression baseline remain supported as legacy demonstrations.
See [ROADMAP.md](ROADMAP.md), [the accepted ADRs](docs/adr/README.md), and
[CONTRIBUTING.md](CONTRIBUTING.md) for scope and evidence rules.

## Wave lab v1.0 quick start

Forward scalar optics uses only the NumPy reference path:

```bash
uv sync --dev
uv run pytest
uv run wave-lab propagate examples/gaussian_lens.toml --out runs/gaussian-lens
```

Differentiable optimization and the full robust flagship use the optional JAX
CPU path in 64-bit mode:

```bash
uv sync --dev --extra jax
JAX_ENABLE_X64=1 uv run --extra jax wave-lab optimize examples/inverse_design_smoke.toml --out runs/inverse-design-smoke
JAX_ENABLE_X64=1 uv run --extra jax wave-lab flagship examples/robust_flagship.toml --out runs/robust-flagship
JAX_ENABLE_X64=1 uv run --extra jax pytest
```

On PowerShell, set `$env:JAX_ENABLE_X64 = "1"` before running the JAX commands.
GPU availability is optional and is not required for correctness or CI.

### Product boundary

v1.0 supports monochromatic coherent scalar fields, Gaussian/Hermite-Gaussian
and arbitrary complex inputs, lenses/apertures/phase masks, Fresnel and
band-limited angular-spectrum propagation, normalized optical objectives,
bounded/smoothed/TV/quantized phase constraints, restartable Adam optimization,
robust ensembles, held-out evaluation, and provenance-rich artifacts.

It does not claim polarization or vector-Maxwell accuracy, nonlinear/broadband
or partially coherent propagation, high-index integrated-photonics validity,
manufacturing-ready layout, or mask export.

### Flagship evidence

The four required comparisons use a common 1064 nm grid, aperture, 12 cm total
length, Gaussian source, HG10 target, constraints, metrics, and seed:

| Baseline | Nominal overlap | Held-out mean | Held-out worst | Worst efficiency |
| --- | ---: | ---: | ---: | ---: |
| No masks | ~0 | ~0 | ~0 | ~0 |
| One optimized mask | 0.608300 | 0.583858 | 0.456924 | 0.456871 |
| Three nominal masks | 0.936511 | 0.899484 | 0.762172 | 0.762128 |
| Three robust masks | 0.924105 | 0.887806 | 0.781349 | 0.781315 |

Optimization and held-out sets are disjoint and both cover wavelength,
alignment, phase-depth, and plane-spacing perturbations. Robust training gains
`0.019178` in held-out worst overlap over nominal three-mask training, while
slightly reducing nominal and held-out mean performance; that tradeoff is part
of the result, not hidden.

The independently implemented NumPy BLAS evaluation reports maximum absolute
Fresnel-transfer loss `4.45e-7`, minimum retained sampled spectral power `1.0`,
grid finest-pair change `0.00643`, padding change `2.11e-6`, aperture change
`1.88e-6`, substep change `0`, and padded boundary power `3.27e-8`. All
predeclared gates pass. The trusted regime is limited to the committed
monochromatic coherent scalar free-space sampling, window, aperture, and angular
content. Optional reduced-scale Maxwell/FDTD was not performed because the
millimetre-scale system is outside practical CPU CI scope; ADR 0004 makes
independent scalar-model transfer the release gate.

Evidence and regeneration details:

- [M1 propagation validation](docs/optics-model-validation.md) and
  [`optics_m1.json`](benchmarks/reference/optics_m1.json)
- [M2 objectives/gradients](docs/objectives-and-gradients.md) and
  [`jax_m2.json`](benchmarks/reference/jax_m2.json)
- [M3 restartable optimizer](docs/inverse-design-runner.md)
- [M4 flagship protocol](docs/flagship-experiment.md) and
  [`flagship_m4.json`](benchmarks/reference/flagship_m4.json)
- [M4 model transfer](docs/flagship-model-transfer.md) and
  [`flagship_validation_m4.json`](benchmarks/reference/flagship_validation_m4.json)

## Interactive dashboard

```bash
uv sync --dev --extra ui
uv run --extra ui streamlit run apps/streamlit_dashboard.py
```

The dashboard leads with NumPy forward optics, constrained inverse design, and
the committed validated flagship summary. Long optimization/flagship actions
are loaded lazily and require the JAX extra. The “Legacy Quantum” workspace
retains free-packet, barrier/Zeno, double-slit, sweep, and quantum-validation
controls. Dashboard outputs remain under ignored `runs/dashboard/` and
`reports/dashboard/` directories.

## Preserved quantum dynamics workflows

This repository now includes a maintained Python package for reproducible
2D quantum dynamics experiments. The original C++/Python/GPU scripts remain
in place as legacy reference material from the Scientific Computing Capstone
project, while the new `quantum-dynamics-lab` package is the recommended path
for new runs.

### Legacy quantum CLI

```bash
uv sync --dev
uv run pytest
uv run quantum-lab run examples/free_packet.toml --out runs/free_packet
uv run quantum-lab run examples/barrier_zeno.toml --out runs/barrier_zeno
uv run quantum-lab compare examples/zeno_sweep.toml --out reports/zeno_sweep
uv run quantum-lab compare examples/zeno_research_sweep.toml --out reports/zeno_research_sweep
uv run quantum-lab run examples/double_slit.toml --out runs/double_slit
uv run quantum-lab narrative --zeno reports/zeno_research_sweep --double-slit runs/double_slit/run.npz --out reports/research_narrative
uv run quantum-lab validate examples/validation_suite.toml --out reports/validation_suite
uv run quantum-lab render runs/free_packet/run.npz --out reports/free_packet
```

### Legacy dashboard interpretation

Install the optional UI dependency and launch the local Streamlit research
dashboard:

```bash
uv sync --dev --extra ui
uv run --extra ui streamlit run apps/streamlit_dashboard.py
```

The legacy workspace provides parameter controls for free-packet, barrier/Zeno,
and double-slit experiments, a saved-frame viewer for probability density and
phase, norm and probability diagnostics, and access to generated artifacts.
Double-slit runs compare coherent and which-path densities,
screen profiles, and a global min/max screen-profile contrast. That contrast is
a visualization and regression diagnostic; it has not been independently
validated as a physical interference-visibility measurement. The Zeno Sweep
workspace runs a bounded grid of barrier heights and measurement intervals,
including an unmeasured baseline, and generates the same CSV, heatmap, and HTML
report used by the CLI comparison workflow. The validation workspace runs the
same numerical checks as `quantum-lab validate` and displays their pass/fail
results and convergence plots. Dashboard output is saved under
`runs/dashboard/` and `reports/dashboard/`, which remain ignored by Git.

The new solver uses a NumPy Strang split-step Fourier method with periodic
FFT boundaries, normalized Gaussian wave packets, configurable potentials, and
static report generation. Run artifacts are saved as compressed `.npz` files
with companion `metrics.json`, `summary.png`, `potential.png`, and, when
enabled and FFmpeg is available, `probability.mp4`.

Supported v1 experiments:

- free packet propagation
- barrier tunneling with repeated Zeno-style projection measurements
- coherent versus which-path double-slit propagation

The maintained implementation lives in `src/quantum_dynamics_lab/`; example
experiment configs live in `examples/`; tests live in `tests/`.

For Zeno barrier experiments, the maintained workflow now supports optional
absorbing edge boundaries and explicit detector placement. Repeated
measurements are modeled as no-click projections: probability in the detector
region is accumulated as detector-click probability, removed from the
conditional wavefunction, and the no-click branch is renormalized while
survival weight tracks the unconditional probability. Generated comparison
reports include `comparison.csv`, `comparison.png`, `zeno_transmission_heatmap.png`,
and `report.md`.

This no-click projection is an algorithmic demonstration, not an independently
validated open-system detector model. Its transmission trend is retained as a
quantum v0.1 regression target; by itself it is not research evidence for the
Quantum Zeno effect.

Optional FFT backends can be selected in `[solver]` with `backend = "numpy"`,
`"pyfftw"`, `"cupy"`, or `"auto"`. The default `auto` uses pyFFTW when it is
installed and otherwise falls back to NumPy; CuPy and pyFFTW remain optional
dependencies. Install them with `uv sync --extra fftw` or
`uv sync --extra cuda12` when the local machine supports those runtimes.

The validation workflow writes `validation_metrics.json`,
`validation_report.md`, `validation_report.html`, and `validation_summary.png`.
It checks norm conservation, free Gaussian dispersion against the analytic
width, qualitative barrier and no-click transmission trends, and optional
backend parity. The norm and dispersion cases are numerical/analytic checks;
the barrier and Zeno trend checks freeze current behavior without independently
validating its physical interpretation. The separately reported double-slit
contrast has the same regression-only status.

## Legacy capstone implementation

The material below describes the original C++/FFTW capstone implementation. It
is retained as historical reference and is not the recommended interface for
new runs.

The mechanics of quantum free particles can be hard to study tangibly. With scientific computing, we are able to build a numerically stable simulation environment to test out models.

This program numerically evolves the linear Schrodinger equation for a complex
scalar wavefunction sampled on a finite, periodic grid.

The original motivating question was whether repeated projection operations in
a grid simulation could illustrate Zeno-like suppression of barrier
transmission. The simulation is a demonstration model and does not establish
the behavior of a physical measurement apparatus.

The static potential and wavefunction are sampled on the same two-dimensional
Cartesian grid.

Operator splitting applies the position-space potential phase and the
momentum-space kinetic phase separately. The Hamiltonian implemented here is
linear; there is no nonlinear Hamiltonian term being removed.

## About 2D Quantum Free Particle

The maintained Python package uses a second-order Strang split-step Fourier
method: a half potential phase, a full spectral kinetic phase, and a second half
potential phase. The legacy C++ loop shown below uses a first-order
split-operator Fourier update with full potential and kinetic phases. Neither
implementation is Crank-Nicolson, and the spatial Laplacian is handled
spectrally rather than by finite differences.

Each legacy time step evolves the wavefunction in the position basis and then
uses a Fourier transform to represent it in the momentum basis.

In the momentum basis, the FFTW implementation applies the kinetic-energy phase.

The wavefunction is finally inverse transformed to position space before the
next time step.

<p align="center">
<img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/visualization/larger-output-quantum.gif"/>
</p>

### Wave function in one dimension X

We have a representation free Schrodinger equation:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?i%5Chbar%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%5CPsi%28t%29%3D%20%5Chat%7BH%7D%5CPsi%28t%29"/>
</p>

where:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28x%29%20%3D%20%5Cpsi_%7B%5Cmathbb%7BR%7D%7D%20%28x%29%20&plus;%20i%20%5Cast%20%5Cpsi_%7B%5Cmathbb%7BI%7D%7D%28x%29"/>
</p>

we take the separable Hamiltonian:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?H%3D%5Cfrac%7Bp%5E2%7D%7B2m%7D&plus;V"/>
</p>

so that:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?i%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%20%5CPsi%28t%29%20%3D%20%5Cleft%28%20%5Cfrac%7B%20%5Chat%7Bp%7D%5E2%20%7D%7B2m%7D&plus;%20V%28x%29%20%5Cright%29%20%5CPsi%28t%29"/>
</p>

so that we have our position operator defined by:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Chat%7Bx%7D%7E%3D%7Ei%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20p%7D%7E%3F"/>
</p>

so that we have our momentum operator defined by:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Chat%7Bp%7D%7E%3D%7E-i%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%7E."/>
</p>

The split-operator method alternates between bases so the potential operator is
pointwise in position space and the kinetic operator is pointwise in momentum
space. FFTW supplies the forward and inverse discrete Fourier transforms.

<p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/visualization/figures/img1.png" width=640 height=480 />
</p>



<p align="center">
<img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/media/sample-quantum-particle.gif"/>
</p>
_____________
## DEPENDENCIES:

___________________________________
LINUX (tested on 16.04 Ubuntu)

### OPTION 1 LINUX: APT-GET PACKAGE MANAGER

in terminal type:

     sudo apt-get install ffmpeg
___________________________________
MAC OS X

### OPTION 1 MAC: HOMEBREW PACKAGE MANAGER

in terminal type:

     brew install ffmpeg

### OPTION 2:
Download latest fftw stable release
http://www.fftw.org/download.html

Once you extract the folder from the ".tar" or ".tar.gz" file navigate to into the "fftw-YOUR_VERSION_NUMBER" directory

Once in the directory
in the terminal enter the following three commands:

    ./configure

    make

    make install

debug:
if error you may not have specified "sudo" privelages
___________________________________________________________________

## COMPILING THE C++ SIMULATION CODE

You will then need to specify the linking flags to compile with the
fftw3 library. On Unix systems, link with "-lfftw3 -lm" like in the example below:

ENTER THIS IN THE MAIN DIRECTORY "2D-Quantum-Free-Particle/":

```bash
g++ 2Dparticle.cpp -o 2Dparticle -lfftw3 -lm
```

RUN THE COMPILED SIMULATION

```c++
./2Dparticle
```

OUTPUT SENT TO DIRECTORY "slices"
Created output will be saved as a ".dat" file in this directory: "2D-Quantum-Free-Particle/slices"

 <p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/media/particle_2D_1.gif"/>
 </p>

## CUDA GPU VERSION

The repository also provides a GPU implementation using CUDA in
`gpu/2D_particle_cuda.cu`. It requires the CUDA toolkit and the cuFFT
library. Compile it with:

```bash
nvcc gpu/2D_particle_cuda.cu -lcufft -o gpu/2DparticleCUDA
```

Running this executable produces text output of the final probability at
the center of the grid.

## Modifying the intial parameters

### Simulation Structure
```c++
     N  - single side length of our N by N array using a standard indexing
     t  - current program time
     dt - time interval between program time steps
     t0 - intial program time
     tf - final program time ( simulation ends once t = tf
     L  - Size of our 2D simulation space in program units
```
### Initial Wave Function Conditions

### Data schema:
Wave function is stored on in multi-dimensional arrays where phi and chi represent the position and momentum space representations of our wave-function.

```c++
     phi[2][2][N][N]  this stores the wavefunction in position space
     phi - position space representation of the wavefunction

     chi[2][2][N][N]  this stores the wavefunction in fourier space
     chi - momentum space representation of the wavefunction

     first row: [2] - Used for before and after partial differential equation solving step
     second row: [2] - Real and Imaginary parts of each wave function
     third row: [N] - X dimension of our 2D simulation space
     fourth row: [N] - Y dimension of our 2D simulation space
```

  ### Variables for keeping track of total energy stored in the simualtion.

```c++
      realsum += psi[0][RE][index];
      complexsum += psi[0][IM][index];
      probabilitysum += psi[0][RE][index] * psi[0][RE][index] + psi[0][IM][index] * psi[0][IM][index];
```        
 ### Gaussian Distribution in 2D

```c++
        sigma = ~0.7
        y = (j*dx) - (L/2.0);
        x = (i*dx) - (L/2.0);
        kx = 1.0*3.141592/L;
        ky = 10.0*3.141592/L;
        A  = AMPLITUDE;
```

  Sets gaussian real part of psi:

```c++
        phi[0][RE][j][i] = A*cos(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));
```

  Sets gaussian for the imaginary parts of psi:

```c++
  phi[0][IM][j][i] = A*sin(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));
```

### Setting Up Fourier Plans
Create both backward and forward fftw plans to be used to time evolve our wavefunction

```c++
//Forward DFT
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
plan = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//Backward DFT
in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
plan2 = fftw_plan_dft_2d(N, N, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);
```
## Time Evolution
time evolution begins at t0 and iterates by time-step dt until current time reaches final time tf
```c++   
    while (t<tf)
```

### Position Space Update
Update wavefunction in position space

```c++    
    {
        for (int j = 0; j < N; j++){
            y = (j*dx - (L/2));
            for (int i = 0; i < N; i++){//update phase of position space wavefunction
                x = (i*dx - (L/2));
                psi[1][RE][j][i] = psi[0][RE][j][i] * cos(potential(x) * dt) + psi[0][IM][j][i]*sin(potential(x) * dt);
                psi[1][IM][j][i] = psi[0][IM][j][i] * cos(potential(x) * dt) - psi[0][RE][j][i]*sin(potential(x) * dt);
            }
        }
```

### Forward Fourier Transform
Update the momentum space representation of the wavefunction
Position Space to momentum space

Load the FFtw arrays for Real and Complex parts of the wavefunction
```c++    
        for (int j = 0; j < N-1; j++){
            for (int i = 0; i < N; i++){  //load our FFTw array
                in[i+j*N][0] = psi[1][RE][j][i];
                in[i+j*N][1] = psi[1][IM][j][i];
            }
        }
```
Execute the plan:
This loop puts the transformed array in DFT output
```c++
        fftw_execute(plan);//transform now stored in out in FFTw format
```       
Unload the output of our fourier transform to chi array where momentum space wavefunction is defined

```c++
            for (int j = 0; j < N-1; j++){
                for (int i = 0; i < N; i++){
                    chi[0][RE][j][i] = out[i+j*N][0];
                    chi[0][IM][j][i] = out[i+j*N][1];
                }
            }
```

### Momentum Space Update
Update the momentum space representation of our wavefunction
```c++    
        for (int j = 0; j < N-1; j++)
        {
            py = ((2*3.145926535)/L) * (( (j + (N/2)) % N) - N/2);
            for (int i = 0; i < N; i++)//here we update the phases in momentum space
            {
                px = ((2*3.145926535)/L) * (( (i + (N/2)) % N) - N/2);
                chi[1][RE][j][i] = chi[0][IM][j][i]*sin((dt*(px*px+py*py))/2) + chi[0][RE][j][i]*cos((dt*(px*px+py*py))/2);
                chi[1][IM][j][i] = chi[0][IM][j][i]*cos((dt*(px*px+py*py))/2) - chi[0][RE][j][i]*sin((dt*(px*px+py*py))/2);
            }
```

### Backwards Fourier Transform
Update our position space representation of our wavefunction

Load the FFtw arrays for Real and Complex parts of the wavefunction
```c++
        for (int j = 0; j < N-1; j++){
            for (int i = 0; i < N; i++){
                in2[i+j*N][0] = chi[1][RE][j][i];
                in2[i+j*N][1] = chi[1][IM][j][i];
            }
        }
```
Excute Bacwards Plan
```c++
        fftw_execute(plan2);
```

Update the position space representation of the wavefunction
```c++
            for (int j = 0; j < N-1; j++){
                for (int i = 0; i < N; i++){
                    //this loop puts the transformed array in DFT output
                    psi[0][RE][j][i] = out2[i+j*N][0];
                    psi[0][IM][j][i] = out2[i+j*N][1];
                }
            }
```

#### Normalize
this loop accounts for unnormalized DFT after forward and backward transforms
```c++
            for (int j = 0; j < N; j++){
                for (int i = 0; i < N; i++){
                    psi[0][RE][j][i] = psi[0][RE][j][i]/(N*N);
                    psi[0][IM][j][i] = psi[0][IM][j][i]/(N*N);
                }
            }
```

### Update & Repeat
Update our time and generation interators
```c++
        t += dt;
        num++;
        printf("*** program time: %lf \n",t);
```

### Output
```c++
        if (num % div == 0)
        {
            outputfield(slicenum);
            //outputenergy(slicenum);
            slicenum++;
        }
```
where:
Output function save file in ".dat" format

```c++
void outputfield(int first)//outputs the field values
{
    static FILE *slicefield;
    static char name[500];
    sprintf(name,"./slices/slices_fields_%d.dat", first);
    slicefield=fopen(name,"w");
    double psiprob[N][N];

    for (int j = 0; j < N; j++)
    {
        for ( int i = 0 ; i < N; i++)
        {
            psiprob[j][i] = psi[0][RE][j][i] * psi[0][RE][j][i] + psi[0][IM][j][i] * psi[0][IM][j][i];
        }
    }
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N; i++)
      {
          if (i%desample==0){
              fprintf(slicefield,"%d  %d  %lf  %lf  %lf", i , j, psi[0][RE][j][i], psi[0][IM][j][i], psiprob[j][i]);
              fprintf(slicefield,"\n");
          }
        }
     }
    fclose(slicefield);
}
```

This repository has a wiki page!
https://github.com/mauckc/2D-Quantum-Free-Particle/wiki

and a python visualization sub-repository!
https://github.com/mauckc/2D-Quantum-Free-Particle/tree/master/visualization

- Ross Mauck
