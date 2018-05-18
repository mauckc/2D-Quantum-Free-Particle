"# 2D-Quantum-Free-Particle" 
This repository has a wiki page!
https://github.com/mauckc/2D-Quantum-Free-Particle/wiki

and a python visualization sub-repository!
https://github.com/mauckc/2D-Quantum-Free-Particle/tree/master/visualization

<p align="center">
  <img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/media/particle_2D_1.gif"/>
</p>

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
     
     i*(dpsi/dx) = - (1/2)*( dpsi/dx)^2 + U(x)* psi(x)
     where psi(x)= psireal(x) + i*psiimag(x)
     
 This version implements a "split-step" Crank-Nicolson method
 We evolve our wave function in the position basis.
 Then we fourier transform the wavefunction to evolve it in the momentum basis
 
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
