# 2D-Quantum-Free-Particle

The mechanics of quantum free particles can be hard to study tangibly. With scientific computing, we are able to build a numerically stable simulation environment to test out models.

This program numerically integrates the Schrodinger equation on finite complex scalar fields for simulating interactions of quantum particles under varied observation.

The goal: Can the Quantum-Zeno Effect a.k.a. the watch dog effect be modeled on an N-dimensional finite difference point lattice. The effects of repeated observation on a quantum wavefunction tend to resist the effects quantum tunneling.

This simulation is built with a 2D static potential field that the wave-function shares the same 2d xy space points.

The optimizations implemented with our work are geared toward removing the need to solve the non-linear term of the hamiltonian in traditional methods.  This frees the computational demand balance with accuracy (non-diverging data)

## About 2D Quantum Free Particle

This version implements a second-order in time finite difference method known as the "split-step" Crank-Nicolson method. By calculating energy states using the hamiltonian in both position and momentum space, this program is able to achieve numerically stable integration, which is necessary for finite difference methods.

Each time-iteration, the program evolves the wave function in the position basis. Then we apply a Fourier transform to the wave function representation in position space to calculate the wave function representation in momentum space.

Once the momentum space representation has been solved by FFTw we can evolve the non-linear term of the Hamiltonian in the momentum basis/phase-space.

The wavefunction in momentum space is finally reverse Fourier transformed back into position space in order to repeat this integration scheme at the next time step t + dt

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

 This version implements a "split-step" Crank-Nicolson method
 We evolve our wave function in the position basis.
 Then we fourier transform the wavefunction to evolve it in the momentum basis

 Instead of integrating the following the standard position space representation of the equation:

 <p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cbg_black%20i%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%5E%7B2%7D%20&plus;%20U%28x%29%5Cpsi%28x%29"/>
</p>

This method is able to achieve a much higher level of accuracy by spliting our integration into two parts and relying on the flexibility of the FFTw algorithm to handle our computational complexity.

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
