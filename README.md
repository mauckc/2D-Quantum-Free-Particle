"# 2D-Quantum-Free-Particle" 

Dependencies:

Create the output directory "2D-Quantum-Free-Particle/slices"

Download latest fftw stable release
http://www.fftw.org/download.html

Once you extract the folder from the .tar file navigate to into the "fftw-YOUR_VERSION_NUMBER" directory

In command prompt enter:

./configure

make

make install

if error you may not have specified "sudo" privelages

You will then need to specify the linking flags to compile with the
fftw3 library. On Unix systems, link with "-lfftw3 -lm".

g++ 2Dparticle.cpp -o 2Dparticle -lfftw3 -lm

./2Dparticle

output located in the slices directory
