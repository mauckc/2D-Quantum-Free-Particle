#
# Python Implementation of 2D-Quantum-Free-Particle schrodinger solver split step method:
#

from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
#import seaborn as sns
import pyfftw
import os
import subprocess
import string

directory = "outfields"

if not os.path.exists(directory):
    os.makedirs(directory)

def calcsum(field):
    total = 0.0
    for i in range(N):
        for j in range(N):
            total += field[i][j]
    return total

def calcavg(field):
    total = 0.0
    for i in range(N):
        for j in range(N):
            total += field[i][j]
    return total/(N*N)

def sum1(input):
    return sum(map(sum, input))

def max1(input):
    return max(map(max, input))

def xnormed(x):
    return  x / x.max()

def calcPsiProbability(psi):
    psiprob = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            psiprob[j,i] = psi[0][0][j][i] * psi[0][0][j][i] + psi[0][1][j][i] * psi[0][1][j][i]
    return psiprob

def checksum(sumvar, n):
    if sumvar is n:
        checkstate = True
    else:
        checkstate = False

    return checkstate

def potential(x):
    U = 10.0
    return U * ( 0.5 - math.pow(math.cos(6*3.141592*(0.65*x-L/2)/L)/2,2.0));

def yxpotential(y,x):
    U = 10.0
    return U * (((0.5 - math.pow(math.cos(6*3.141592*(0.65*x-L/2)/L)/2,2.0)) + (0.5 - math.pow(math.cos(6*3.141592*(0.65*y-L/2)/L)/2,2.0)))/2) + x + y; # + x adds gradient in x direction

if __name__ == '__main__':
    #x = int(raw_input())
    #y = int(raw_input())
    #z = int(raw_input())
    #n = int(raw_input())

    # define initial parameters
    N = 256 # number of evenly spaced points
    L = 80.0 # Length of the box ( box is simulation space in world coordinates )
    dt = 0.005 # Time-step
    t0 = 0.0 # Initial time
    tf = 2.0 # final time
    dx = L/N # distance between points in program dimensions
    t = 0.0
    desample = 8 # desampling factor for output after calculation
    # amplitude for initial potential settings
    amplitude = 1.0
    A = amplitude
    # fill in the rest of the initial settings
    # sigma of gaussian
    sigma = 5.73
    num = 0 # for outfield iterator
    real = 0 # real and imaginary parts of our wave-function
    imag = 1 # for accessing the array parts with more clarity

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

    # Setup initial gaussian
    for i in range(N):
        y = (i * dx) - (L/3.8) # or 3.8
        for j in range(N):
            x = (j*dx) - (L/3.8) # or 3.8
            kx = ((0.7 * math.pi) / L)
            ky = ((1.0 * math.pi) / L)
            A = amplitude

            # set the gaussian fields in position space
            psi[0][real][i][j] = A * math.cos(kx*x+ky*y) * math.exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)))
            psi[0][imag][i][j] = A * math.sin(kx*x+ky*y) * math.exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)))

            phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]


    print("Initial Conditions Set..")

    # begin time evolution
    plt.imshow(phi,cmap='hot',interpolation='nearest')
    plt.savefig('outfields/field%04d.png' % num)
    plt.show()
    num += 1
    outnum = 1
    # alternative ploting
    #ax = sns.heatmap(phi, linewidth=0.5)

    # Start time simulation loop
    while(t<tf):

        # update field in position space with respect to potential for real
        #  and imaginary parts
        for i in range(N):
            y = (i*dx - (L/2))
            for j in range(N):
                x = (j*dx - (L/2))
                psi[1][real][i][j] = psi[0][real][i][j] * math.cos(yxpotential(y,x) * dt) + psi[0][imag][i][j] * math.sin(yxpotential(y,x)*dt)
                psi[1][imag][i][j] = psi[0][imag][i][j] * math.cos(yxpotential(y,x) * dt) + psi[0][real][i][j] * math.sin(yxpotential(y,x)*dt)
        # load the input array for the transform
        for i in range(N):
            for j in range(N):
                forwardin[i,j] = psi[1][real][i][j] + 1j*psi[1][imag][i][j]
        # perform the transform
        forwardout = pyfftw.interfaces.numpy_fft.fft(forwardin)

        #print(forwardout)
        #print(np.shape(forwardout))
        for i in range(N):
            for j in range(N):
                chi[0][real][i][j] = forwardout[i,j].real
                chi[0][imag][i][j] = forwardout[i,j].imag

        # Nomalize before updating in momentum space
        chi[0][real] = xnormed(chi[0][real])
        chi[0][imag] = xnormed(chi[0][imag])
        # update field in momentum space
        for i in range(N):
            py = ((2*math.pi)/L) * (((i + (N/2)) / N) - N/2)
            for j in range(N):
                px = ((2*math.pi)/L) * (((j + (N/2)) / N) - N/2)
                chi[1][real][i][j] = chi[0][imag][i][j] * math.sin((dt*(px*px+py*py))/2) + chi[0][real][i][j] * math.cos((dt*(px*px+py*py))/2)
                chi[1][imag][i][j] = chi[0][imag][i][j] * math.cos((dt*(px*px+py*py))/2) - chi[0][real][i][j] * math.sin((dt*(px*px+py*py))/2)

        # This could be prettified using list comprehensions?
        for i in range(N):
            for j in range(N):
                backwardin[i,j] = chi[1][real][i][j] + 1j*chi[1][imag][i][j]

        # switch back to position spaces
        backwardout = pyfftw.interfaces.numpy_fft.fft(backwardin)
        for i in range(N):
            for j in range(N):
                psi[0][real][i][j] = backwardout[i,j].real
                psi[0][imag][i][j] = backwardout[i,j].imag

        psi[0][real] = xnormed(psi[0][real])
        psi[0][imag] = xnormed(psi[0][imag])

        # for i in range(N):
        #     for j in range(N):
        #         phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]

        phi = calcPsiProbability(psi)
        # normalize phi
        #phi = xnormed(phi)

        print("Phi sum:"+str(calcsum(phi)))
        print("Phi sum1:"+str(sum1(phi)))
        print("Phi max:"+str(max1(phi)))
        print("Phi avg:"+str(calcavg(phi))+"\n")
        if num % desample is 0:
            #phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
            plt.imshow(phi,cmap='hot',interpolation='nearest')
            plt.savefig('outfields/field%04d.png' % outnum)
            #plt.show()

            print("printing field number "+str(num)+" output number "+str(outnum))
            print( "at time: "+str(t))
            outnum += 1
        num += 1

        t += dt
    # when simulation has reached time final
    #plt.imshow(phi,cmap='hot',interpolation='nearest')
    #plt.show()
    outdt = str(dt).translate(None, string.punctuation)
    returncode = subprocess.call(['ffmpeg','-r','30','-f','image2','-s','640x480','-i','outfields/field%'+'04d.png','-vcodec','libx264','-crf','25','-pix_fmt','yuv420p','fields_out_N%d_L%d_dt%s_dsmp%d_tf%d_sigma%d.mp4' % (N,math.floor(L),outdt,desample, math.floor(tf), math.floor(sigma))])
    print('returncode:', returncode)
    # call(["./argpngs2mp4.sh %s %s %s %s " % (str(N),str(round(L,1)),str(round(dt,3)),str(desample))])
