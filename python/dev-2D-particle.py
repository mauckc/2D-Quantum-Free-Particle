#
# Python Implementation of 2D-Quantum-Free-Particle schrodinger solver split step method:
#

from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import pyfftw

def calcsum(field):
    total = 0.0
    for i in range(N):
        for j in range(N):
            total += field[i][j]
    return total

def calcPsiProbability(phi, chi):
    for i in range(N):
        for j in range(N):
            psiprob[i,j] = psi[0][real][j][i] * psi[0][real][j][i] + psi[0][imag][j][i] * psi[0][imag][j][i]
    return psiprob

def outputfield(outfieldnumber):
    # open file
    # output each field like c++ i.e.: fprintf(slicefield,"%d %d %lf %lf", i, j, psi[0][real][i][j],psi[0][imag][i][j], psiprob[i][j])
    return

def ddfield(field, dfield):

    return # second derivative of the field

def sumxyz(x, y, z):
    sumvar = x + y + z
    return sumvar

def checksum(sumvar, n):
    if sumvar is n:
        checkstate = True
    else:
        checkstate = False

    return checkstate

L = 5.0 # Length of the box ( box is simulation space in world coordinates )

def potential(x):
    U = 50.0
    return U * ( 1.0 - math.pow(math.cos(6*3.141592*(0.65*x-L/2)/L)/2,2.0));

if __name__ == '__main__':
    #x = int(raw_input())
    #y = int(raw_input())
    #z = int(raw_input())
    #n = int(raw_input())

    # Create array of the input point to print them out
    #ar = [x,y,z]
    #print(ar)
    # define initial parameters
    N = 256 # number of evenly spaced points
    L = 20.0 # Length of the box ( box is simulation space in world coordinates )
    dt = 0.01 # Time-step
    t0 = 0.0 # Initial time
    tf = 5.0 # final time
    dx = L/N # distance between points in program dimensions

    num = 0 # for outfield iterator

    #psi[2][2][N][N]  # wavefunction array in position space
    #chi[2][2][N][N] # wavefunction in momentum space

    div = 1 # scale factor for reducing frequency of the time-step
    desample = 4 # desampling factor for output after calculation

    # amplitude for initial potential settings
    amplitude = 1.0
    A = 1.
    # fill in the rest of the initial settings
    # sigma of gaussian
    sigma = 0.73

    realsum = 0.0
    complexsum = 0.0
    probabilitysum = 0.0
    slicenum = 0 # count number of slices

    #x = 0.0
    #y = 0.0
    #p = 0.0
    #px = 0.0
    #py = 0.0
    #kx = 0.0
    #ky = 0.0

    num = 0 # used for output iterator
    real = 0 # real and imaginary parts of our wave-function
    imag = 1 # for accessing the array parts with more clarity


    t = 0.0
    # initialize arrays
    psi = np.zeros((2,2,N,N)) # wavefunction in  position space
    chi = np.zeros((2,2,N,N)) # wavefunction in fourier space
    phi = np.zeros((N,N)) # probability of wavefunctions

    print("fields structure allocated in position space -> dimensions:")
    print(np.shape(phi))
    print("fields structure allocated in position space -> dimensions:")
    print(np.shape(chi))
    print("fields structure allocated in probability space -> dimensions:")
    print(np.shape(psi))


    # print initial settings and setup information
    # initial fftw plans if implemented later see C++ implemenation for structure
    # Forward DFT
    forwardin = pyfftw.empty_aligned((N,N), dtype='complex128')
    forwardout = pyfftw.empty_aligned((N,N), dtype='complex128')

    # Backward DFT
    backwardin = pyfftw.empty_aligned((N,N), dtype='complex128')
    backwardout = pyfftw.empty_aligned((N,N), dtype='complex128')


    for i in range(N):
        y = (i * dx) - (L/3.8)
        for j in range(N):
            x = (i*dx) - (L/3.8)
            kx = ((20.0 * math.pi) / L)
            ky = ((4.0 * math.pi) / L)
            A = amplitude

            # set the gaussian fields in position space
            psi[0][real][i][j] = A * math.cos(kx*x+ky*y) * math.exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)))
            psi[0][imag][i][j] = A * math.sin(kx*x+ky*y) * math.exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)))

            phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
    print("Initial Conditions Set..")

    #print(phi)

    # begin time evolution
    plt.imshow(phi,cmap='hot',interpolation='nearest')
    plt.show()

    # alternative ploting
    ax = sns.heatmap(phi, linewidth=0.5)
    plt.show()

    while(t<tf):

        # update field in position space with respect to potential for real
        #  and imaginary parts
        for i in range(N):
            y = (i*dx - (L/2))
            for j in range(N):
                x = (j*dx - (L/2))
                psi[1][real][i][j] = psi[0][real][i][j] * math.cos(potential(x) * dt) + psi[0][imag][i][j] * math.sin(potential(x)*dt)
                psi[1][imag][i][j] = psi[0][imag][i][j] * math.cos(potential(x) * dt) + psi[0][real][i][j] * math.sin(potential(x)*dt)
        # load the input array for the transform
        for i in range(N):
            for j in range(N):
                forwardin[i,j] = psi[1][real][i][j] + 1j*psi[1][imag][i][j]
        # perform the transform
        forwardout = pyfftw.interfaces.numpy_fft.fft(forwardin)

        #print(forwardout)
        print(np.shape(forwardout))
        for i in range(N):
            for j in range(N):
                chi[0][real][i][j] = forwardout[i,j].real
                chi[0][imag][i][j] = forwardout[i,j].imag

        # update field in momentum space
        for i in range(N):
            px = ((2*math.pi)/L) * (((i + (N/2)) / N) - N/2)
            for j in range(N):
                py = ((2*math.pi)/L) * (((j + (N/2)) / N) - N/2)
                chi[1][real][i][j] = chi[0][imag][i][j] * math.sin((dt*(px*px+py*py))/2) + chi[0][real][i][j] * math.cos((dt*(px*px+py*py))/2)
                chi[1][imag][i][j] = chi[0][imag][i][j] * math.cos((dt*(px*px+py*py))/2) - chi[0][real][i][j] * math.sin((dt*(px*px+py*py))/2)

        for i in range(N):
            for j in range(N):
                backwardin[i,j] = chi[1][real][i][j] + 1j*chi[1][imag][i][j]

        # switch back to position spaces

        backwardout = pyfftw.interfaces.numpy_fft.fft(backwardin)
        for i in range(N):
            for j in range(N):
                psi[0][real][i][j] = backwardout[i,j].real
                psi[0][imag][i][j] = backwardout[i,j].imag

        # piq
        if num % 10 is 0:
            phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
            plt.imshow(phi,cmap='hot',interpolation='nearest')
            plt.show()
            print("\n printing field number ")
            print( num )
            print( " at time: ")
            print(t)
        num += 1
        t = t + dt
    # when simulation has reached time final
    plt.imshow(phi,cmap='hot',interpolation='nearest')
    plt.show()
