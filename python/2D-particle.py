#
# Python Implementation of 2D-Quantum-Free-Particle schrodinger solver split step method:
#

from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
import pyfftw
import os

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
    if sumvar == n:
        checkstate = True
    else:
        checkstate = False

    return checkstate

def potential(x):
    U = 40.0
    return U * (1.0 - np.power(np.cos(6 * np.pi * (0.65 * x - L / 2) / L) / 2, 2.0))

# def double_slit_potential(x, y, L, barrier_width, slit_width, barrier_potential):
#     # The barrier is in the middle of the domain
#     barrier_center = L / 2
#     slit_center = L / 2
#     # Check if the point is within the barrier region
#     if barrier_center - barrier_width / 2 < x < barrier_center + barrier_width / 2:
#         # Check if the point is within either of the two slits
#         if slit_center - 3 * slit_width / 2 < y < slit_center - slit_width / 2 or \
#            slit_center + slit_width / 2 < y < slit_center + 3 * slit_width / 2:
#             return 0.0
#         else:
#             return barrier_potential
#     else:
#         return 0.0

def double_slit_potential(x, y, L, barrier_width, slit_width, slit_distance, barrier_potential):
    barrier_center = L / 2
    slit_center = L / 2

    if barrier_center - barrier_width / 2 < x < barrier_center + barrier_width / 2:
        half_slit_length = slit_width / 2
        half_slit_distance = slit_distance / 2
        
        if slit_center - half_slit_length - half_slit_distance < y < slit_center - half_slit_length + half_slit_distance or \
           slit_center + half_slit_length - half_slit_distance < y < slit_center + half_slit_length + half_slit_distance:
            return 0.0
        else:
            return barrier_potential
    else:
        return 0.0

def visualize_potential(L, N):
    x_values = np.linspace(-L / 2, L / 2, N)
    y_values = np.linspace(-L / 2, L / 2, N)
    X, Y = np.meshgrid(x_values, y_values)
    Z = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            Z[i, j] = potential(X[i, j])

    plt.figure()
    plt.imshow(Z, cmap='viridis', extent=(-L / 2, L / 2, -L / 2, L / 2), origin='lower')
    plt.colorbar(label='Potential')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Potential Function Visualization')
    plt.show()

def visualize_double_slit_potential(L, N, barrier_width, slit_width, slit_distance, barrier_potential):
    x_values = np.linspace(0, L, N)
    y_values = np.linspace(0, L, N)
    X, Y = np.meshgrid(x_values, y_values)
    Z = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            Z[i, j] = double_slit_potential(X[i, j], Y[i, j], L, barrier_width, slit_width, slit_distance, barrier_potential)

    plt.figure()
    plt.imshow(Z, cmap='viridis', extent=(0, L, 0, L), origin='lower')
    plt.colorbar(label='Potential')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Double-slit Potential Function Visualization')
    plt.show()

if __name__ == '__main__':
    #x = int(raw_input())
    #y = int(raw_input())
    #z = int(raw_input())
    #n = int(raw_input())

    # Create array of the input point to print them out
    #ar = [x,y,z]
    #print(ar)
    # define initial parameters
    N = 128 # number of evenly spaced points
    L = 10.0 # Length of the box ( box is simulation space in world coordinates )
    dt = 0.001 # Time-step
    t0 = 0.0 # Initial time
    tf = 5.0 # final time
    dx = L/N # distance between points in program dimensions

    barrier_width = 0.5
    slit_width = 1.0
    slit_distance = 0.5
    barrier_potential = 200.0


    
    num = 0 # for outfield iterator

    #psi[2][2][N][N]  # wavefunction array in position space
    #chi[2][2][N][N] # wavefunction in momentum space

    div = 1 # scale factor for reducing frequency of the time-step
    desample = 4 # desampling factor for output after calculation

    kx_magnitude = -5.0
    ky_magnitude = 0.0
    # amplitude for initial potential settings
    amplitude = 1.0
    A = amplitude
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
        y = (i * dx) - (L/2)
        for j in range(N):
            x = (j*dx) - (L/3.8)
            kx = ((kx_magnitude * math.pi) / L)
            ky = ((ky_magnitude * math.pi) / L)
            A = amplitude

            # set the gaussian fields in position space
            psi[0][real][i][j] = A * math.cos(kx*x+ky*y) * math.exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)))
            psi[0][imag][i][j] = A * math.sin(kx*x+ky*y) * math.exp(-((x*x)/(4*sigma*sigma) + (y*y)/(4*sigma*sigma)))

            phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
    

    print("Visualizing the potential")
    # visualize_potential(L, N)
    visualize_double_slit_potential(L, N, barrier_width, slit_width, slit_distance, barrier_potential)

    print("Initial Conditions Set..")

    # begin time evolution
    plt.imshow(phi,cmap='hot',interpolation='nearest')
    plt.show()
    if not os.path.exists('outfields'):
        os.mkdir('outfields')
    plt.savefig('outfields/field%04d.png' % num)
    num += 1
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
                psi[1][real][i][j] = psi[0][real][i][j] * math.cos(double_slit_potential(x,y,L,barrier_width,slit_width,slit_distance,barrier_potential) * dt) + psi[0][imag][i][j] * math.sin(double_slit_potential(x,y,L,barrier_width,slit_width,slit_distance,barrier_potential)*dt)
                psi[1][imag][i][j] = psi[0][imag][i][j] * math.cos(double_slit_potential(x,y,L,barrier_width,slit_width,slit_distance,barrier_potential) * dt) + psi[0][real][i][j] * math.sin(double_slit_potential(x,y,L,barrier_width,slit_width,slit_distance,barrier_potential)*dt)
        # load the input array for the transform
        for i in range(N):
            for j in range(N):
                forwardin[i,j] = psi[1][real][i][j] + 1j*psi[1][imag][i][j]
        # perform the transform
        forwardout = pyfftw.interfaces.numpy_fft.fft2(forwardin)

        #print(forwardout)
        print(np.shape(forwardout))
        for i in range(N):
            for j in range(N):
                chi[0][real][i][j] = forwardout[i,j].real
                chi[0][imag][i][j] = forwardout[i,j].imag

        k_values = np.fft.fftfreq(N, d=dx) * 2 * np.pi
        k_values_shifted = np.fft.ifftshift(k_values)
        # update field in momentum space
        for i in range(N):
            py = k_values_shifted[i]
            for j in range(N):
                px = k_values_shifted[j]
                chi[1][real][i][j] = chi[0][imag][i][j] * math.sin((dt*(px*px+py*py))/2) + chi[0][real][i][j] * math.cos((dt*(px*px+py*py))/2)
                chi[1][imag][i][j] = chi[0][imag][i][j] * math.cos((dt*(px*px+py*py))/2) - chi[0][real][i][j] * math.sin((dt*(px*px+py*py))/2)

        # This could be prettified using list comprehensions?
        for i in range(N):
            for j in range(N):
                backwardin[i,j] = chi[1][real][i][j] + 1j*chi[1][imag][i][j]

        # switch back to position space
        backwardout = pyfftw.interfaces.numpy_fft.ifft2(backwardin)

        for i in range(N):
            for j in range(N):
                psi[0][real][i][j] = backwardout[i,j].real
                psi[0][imag][i][j] = backwardout[i,j].imag

        for i in range(N):
            for j in range(N):
                phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]

        print("Phi sum:")
        print(calcsum(phi))
        if num % div == 0:
            #phi[i,j] = psi[0][real][i][j] * psi[0][real][i][j] + psi[0][imag][i][j] * psi[0][imag][i][j]
            plt.imshow(phi,cmap='hot',interpolation='nearest')
            plt.savefig('outfields/field%04d.png' % num)
            plt.pause(.05)
            plt.draw()
            # plt.show()

            print("\nprinting field number ")
            print( num )
            print( "at time: ")
            print(t)
        num += 1
        t += dt
    # when simulation has reached time final
    plt.imshow(phi,cmap='hot',interpolation='nearest')
    plt.show()
