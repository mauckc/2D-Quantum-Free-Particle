# 2D-Quantum-Free-Particle

<p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/media/particle_2D_1.gif"/>
</p>

This program numerically integrates the Schrodinger equation on finite complex scalar fields for simulating interactions of quantum particles under varied observation.

This version implements a second-order in time finite difference method known as the "split-step" Crank-Nicolson method. By calculating energy states using the hamiltonian in both position and momentum space, this program is able to achieve numerically stable integration, which is necessary for finite difference methods.

Each time-iteration, the program evolves the wave function in the position basis. Then we apply a Fourier transform the wave function to evolve the non-linear term of the Hamiltonian in the momentum basis/phase-space. The waveform is then reverse Fourier transformed back into position space in order to repeat this evolution of the waveform.

This repository has a wiki page!
https://github.com/mauckc/2D-Quantum-Free-Particle/wiki

and a python visualization sub-repository!
https://github.com/mauckc/2D-Quantum-Free-Particle/tree/master/visualization

_________________________________________________
HOW TO MODIFY PROGRAM SPACE AND INITIAL CONDITION

 Simulation Structure
 
     N - single side length of our N by N array using a standard indexing
     t - current program time
     dt - time interval between program time steps
     t0 - intial program time
     tf - final program time ( simulation ends once t = tf
     L - Size of our 2D simulation space in program units
 
Initial Wave Function Conditions 
 phi - Real Part of the wavefunction
 chi - Complex Part of the wavefunction

     PARTIAL DIFFERNTIAL EQUATION EXAMPLE FOR WAVEFUNTION IN ONE DIMENSION (X)

  <img src="https://latex.codecogs.com/gif.latex?%5Cbg_black%20i%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%5E%7B2%7D%20&plus;%20U%28x%29%5Cpsi%28x%29"/>
</p>

where:

<img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28x%29%20%3D%20%5Cpsi_%7B%5Cmathbb%7BR%7D%7D%20%28x%29%20&plus;%20i%20%5Cast%20%5Cpsi_%7B%5Cmathbb%7BI%7D%7D%28x%29"/>
</p>
     
 This version implements a "split-step" Crank-Nicolson method
 We evolve our wave function in the position basis.
 Then we fourier transform the wavefunction to evolve it in the momentum basis
 
 <p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/visualization/figures/img1.png" width=640 height=480 />
</p>
 
 data schema:
 
     phi[2][2][N][N]
     chi[2][2][N][N]
 
     first row: [2] - Used for before and after partial differential equation solving step
     second row: [2] - Real and Imaginary parts of each wave function
     third row: [N] - X dimension of our 2D simulation space
     fourth row: [N] - Y dimension of our 2D simulation space
 
  Variables for keeping track of total energy stored in the simualtion.
  
      realsum += psi[0][RE][index];
      complexsum += psi[0][IM][index];
      probabilitysum += psi[0][RE][index] * psi[0][RE][index] + psi[0][IM][index] * psi[0][IM][index];
                  
 Gaussian Distribution in 2D
        
        sigma = ~0.7
        y = (j*dx) - (L/2.0);
        x = (i*dx) - (L/2.0);
        kx = 1.0*3.141592/L;
        ky = 10.0*3.141592/L;
        A  = AMPLITUDE;
        
  Sets gaussian real part of psi
        
        psi[0][RE][j][i] = A*cos(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));
        
  Sets gaussian for the imaginary parts of psi
  
        psi[0][IM][j][i] = A*sin(kx*x+ky*y) * exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)));
        
_____________
DEPENDENCIES:

___________________________________
LINUX (tested on 16.04 Ubuntu)

OPTION 1: APT-GET PACKAGE MANAGER

in terminal type: 

     sudo apt-get install ffmpeg
___________________________________
MAC OS X (tested on 16.04 Ubuntu)

OPTION 1: HOMEBREW PACKAGE MANAGER

in terminal type: 

     brew install ffmpeg

OPTION 2:
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
COMPILING THE C++ SIMULATION CODE

You will then need to specify the linking flags to compile with the
fftw3 library. On Unix systems, link with "-lfftw3 -lm" like in the example below:

ENTER THIS IN THE MAIN DIRECTORY "2D-Quantum-Free-Particle/":

      g++ 2Dparticle.cpp -o 2Dparticle -lfftw3 -lm

RUN THE COMPILED SIMULATION

      ./2Dparticle

OUTPUT SENT TO DIRECTORY "slices"
Created output will be saved as a ".dat" file in this directory: "2D-Quantum-Free-Particle/slices"
