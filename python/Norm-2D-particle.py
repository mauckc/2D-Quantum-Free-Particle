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

def calcPsiProbability(phi, chi):
    for i in range(N):
        for j in range(N):
            psiprob[i,j] = psi[0][real][j][i] * psi[0][real][j][i] + psi[0][imag][j][i] * psi[0][imag][j][i]
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
    return U * (((0.5 - math.pow(math.cos(6*3.141592*(0.65*x-L/2)/L)/2,2.0)) + (0.5 - math.pow(math.cos(6*3.141592*(0.65*y-L/2)/L)/2,2.0)))/2);

if __name__ == '__main__':
    #x = int(raw_input())
    #y = int(raw_input())
    #z = int(raw_input())
    #n = int(raw_input())

    #psi[2][2][N][N]  # wavefunction array in position space
    #chi[2][2][N][N] # wavefunction in momentum space
    # Create array of the input point to print them out
    #ar = [x,y,z]
    #print(ar)
    # define initial parameters
    N = 128 # number of evenly spaced points
    L = 60.0 # Length of the box ( box is simulation space in world coordinates )
    dt = 0.01 # Time-step
    t0 = 0.0 # Initial time
    tf = 3.0 # final time
    dx = L/N # distance between points in program dimensions

    num = 0 # for outfield iterator
    div = 4 # scale factor for reducing frequency of the time-step
    desample = div # desampling factor for output after calculation

    # amplitude for initial potential settings
    amplitude = 2.0
    A = amplitude
    # fill in the rest of the initial settings
    # sigma of gaussian
    sigma = 5.73

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
        y = (i * dx) - (L/3.8) # or 3.8
        for j in range(N):
            x = (j*dx) - (L/3.8) # or 3.8
            kx = ((1000.7 * math.pi) / L)
            ky = ((10.0 * math.pi) / L)
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
        print(np.shape(forwardout))
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

        for i in range(N):
            for j in range(N):
                phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
        # normalize phi
        #phi = xnormed(phi)

        print("Phi sum:"+str(calcsum(phi)))
        print("Phi sum1:"+str(sum1(phi)))
        print("Phi max:"+str(max1(phi)))
        print("Phi avg:"+str(calcavg(phi)))
        if num % desample is 0:
            #phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
            plt.imshow(phi,cmap='hot',interpolation='nearest')
            plt.savefig('outfields/field%04d.png' % outnum)
            #plt.show()

            print("\nprinting field number "+str(num)+" output number "+str(outnum))
            print( "at time: "+str(t))
            outnum += 1
        num += 1

        t += dt
    # when simulation has reached time final
    plt.imshow(phi,cmap='hot',interpolation='nearest')
    plt.show()
    N = 128 # number of evenly spaced points
    L = 60.0 # Length of the box ( box is simulation space in world coordinates ) 60.0
    dt = 0.01 # Time-step 0.01
    t0 = 0.0 # Initial time
    tf = 50.0 # final time 50.0
    dx = L/N # distance between points in program dimensions

    num = 0 # for outfield iterator
    div = 4 # scale factor for reducing frequency of the time-step
    desample = div # desampling factor for output after calculation

    # amplitude for initial potential settings
    amplitude = 1.0 # 2.0
    A = amplitude
    # fill in the rest of the initial settings
    # sigma of gaussian
    sigma = 5.73 # 5.73
    outdt = str(dt).translate(None, string.punctuation)
    returncode = subprocess.call(['ffmpeg','-r','30','-f','image2','-s','640x480','-i','outfields/field%'+'04d.png','-vcodec','libx264','-crf','25','-pix_fmt','yuv420p','fields_out_%d_%d_%s_%d.mp4' % (N,math.floor(L),outdt,desample)])
    print('returncode:', returncode)
    # call(["./argpngs2mp4.sh %s %s %s %s " % (str(N),str(round(L,1)),str(round(dt,3)),str(desample))])
