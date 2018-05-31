# 2D-Quantum-Free-Particle

This program numerically integrates the Schrodinger equation on finite complex scalar fields for simulating interactions of quantum particles under varied observation.

This repository has a wiki page!
https://github.com/mauckc/2D-Quantum-Free-Particle/wiki

and a python visualization sub-repository!
https://github.com/mauckc/2D-Quantum-Free-Particle/tree/master/visualization

<p align="center">
<img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/visualization/larger-output-quantum.gif"/>
</p>

## About 2D Quantum Free Particle

This version implements a second-order in time finite difference method known as the "split-step" Crank-Nicolson method. By calculating energy states using the hamiltonian in both position and momentum space, this program is able to achieve numerically stable integration, which is necessary for finite difference methods.

Each time-iteration, the program evolves the wave function in the position basis. Then we apply a Fourier transform to the wave function representation in position space to calculate the wave function representation in momentum space.

Once the momentum space representation has been solved by FFTw we can evolve the non-linear term of the Hamiltonian in the momentum basis/phase-space. 

The wavefunction in momentum space is finally reverse Fourier transformed back into position space in order to repeat this integration scheme at the next time step t + dt

### PARTIAL DIFFERNTIAL EQUATION EXAMPLE 
for wave function in one dimension X

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cbg_black%20i%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%5E%7B2%7D%20&plus;%20U%28x%29%5Cpsi%28x%29"/>
</p>

where:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28x%29%20%3D%20%5Cpsi_%7B%5Cmathbb%7BR%7D%7D%20%28x%29%20&plus;%20i%20%5Cast%20%5Cpsi_%7B%5Cmathbb%7BI%7D%7D%28x%29"/>
</p>

 We have a representation free Schrodinger equation:
 
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?i%5Chbar%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%5CPsi%28t%29%3D%20%5Chat%7BH%7D%5CPsi%28t%29"/>
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
 

<p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/visualization/figures/img1.png" width=640 height=480 />
</p>

 
## Modifying the intial parameters

### Simulation Structure
 
     N  - single side length of our N by N array using a standard indexing
     t  - current program time
     dt - time interval between program time steps
     t0 - intial program time
     tf - final program time ( simulation ends once t = tf
     L  - Size of our 2D simulation space in program units
     
### Initial Wave Function Conditions 

### Data schema:

     phi[2][2][N][N]  this stores the wavefunction in position space 
     phi - position space representation of the wavefunction

     chi[2][2][N][N]  this stores the wavefunction in fourier space
     chi - momentum space representation of the wavefunction
 
     first row: [2] - Used for before and after partial differential equation solving step
     second row: [2] - Real and Imaginary parts of each wave function
     third row: [N] - X dimension of our 2D simulation space
     fourth row: [N] - Y dimension of our 2D simulation space
 
  ### Variables for keeping track of total energy stored in the simualtion.
  
      realsum += psi[0][RE][index];
      complexsum += psi[0][IM][index];
      probabilitysum += psi[0][RE][index] * psi[0][RE][index] + psi[0][IM][index] * psi[0][IM][index];
                  
 ### Gaussian Distribution in 2D
        
        sigma = ~0.7
        y = (j*dx) - (L/2.0);
        x = (i*dx) - (L/2.0);
        kx = 1.0*3.141592/L;
        ky = 10.0*3.141592/L;
        A  = AMPLITUDE;
        
  Sets gaussian real part of psi:
        
        phi[0][RE][j][i] = A*cos(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));
        
  Sets gaussian for the imaginary parts of psi:
  
        phi[0][IM][j][i] = A*sin(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));
        
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

 <p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/media/particle_2D_1.gif"/>
</p>

## COMPILING THE C++ SIMULATION CODE

You will then need to specify the linking flags to compile with the
fftw3 library. On Unix systems, link with "-lfftw3 -lm" like in the example below:

ENTER THIS IN THE MAIN DIRECTORY "2D-Quantum-Free-Particle/":

      g++ 2Dparticle.cpp -o 2Dparticle -lfftw3 -lm

RUN THE COMPILED SIMULATION

      ./2Dparticle

OUTPUT SENT TO DIRECTORY "slices"
Created output will be saved as a ".dat" file in this directory: "2D-Quantum-Free-Particle/slices"

- Ross Mauck


