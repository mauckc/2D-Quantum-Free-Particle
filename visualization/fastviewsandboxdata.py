# load and make a slice png

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np
from pylab import *
import subprocess
from datetime import datetime

# set starttime values
t1 = datetime.now()

# Desample rate from c++ output
desampled = 4 # 4

# Set dimension of our 2D grid
n = 256/(desampled/2) #256
noutputs = 200
# Boolean to choose to rotate the output graph
rotate = True




print('Setting Dimensions\n n=%d \n' % (n))
print('noutputs=%d \n' % (noutputs))

# Load data files from slices output directory
DataIn = np.loadtxt('slices/slices_fields_0.dat')

# Create figure for plot
fig = plt.figure(figsize=(8,8))

# Populate 1-Ds arrays with length n of floating point values between 0 and 1 for our bounds
x = np.linspace(0., 1., n)
y = np.linspace(0., 1., n)
#The numpy.meshgrid() function generates the samples from an explicit 2D function.
X, Y = np.meshgrid(x, y)
Z = np.zeros((n, n))

# Create array for storing each field dat file data
AllSlicesData = np.zeros((noutputs, (n*n), 5))
Allprob = np.zeros((noutputs, (n*n) ))

# For number of outputs user variable load a successive slice
for k in xrange(noutputs):
    AllSlicesData[k] = np.loadtxt('slices/slices_fields_%d.dat' % (k))
    Allprob[k] = AllSlicesData[k,::,4]
    if k % (4) == 0:
        print('loading slice %d...\n' % (k))

# Interpret All slices data fields
    #indexA = AllSlicesData[-1,::,0]
    #indexB = AllSlicesData[-1,::,1]
    #phi = AllSlicesData[-1,::,2]
    #chi = AllSlicesData[-1,::,3]
    #prob = AllSlicesData[-1,::,4]

#for p in xrange(10):
    #print Allprob[p]

# Interpret single slice data fields
#indexA = DataIn[::,0]
#indexB = DataIn[::,1]
#phi = DataIn[::,2]
#chi = DataIn[::,3]
#prob = DataIn[::,4]

#for k in xrange(noutputs):
    # Must interpret our indexed array and save it to the new data type
    for j in xrange(0, n-1):
        for i in xrange(0, n-1):
            Z[i][j] = Allprob[k,i+j*(n-1)]

    # Normalize to [0,1]
    Z = (Z-Z.min())/(Z.max()-Z.min())

    #colors = cm.viridis(Z)
    #rcount, ccount, _ = colors.shape

    fig = plt.figure()

    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    ax1.plot_wireframe(X, Y, Z, rstride=(desampled/(8/desampled)), cstride=(desampled/(8/desampled)))

    #ax = fig.add_subplot(111, projection='3d')
    #X, Y, Z = axes3d.get_test_data(0.1)

    #ax1 = fig.gca(projection='3d')

    #surf = ax1.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount, facecolors=colors, shade=False)

    ax1.set_xlabel('x - dimension')
    ax1.set_ylabel('y - spatial dimension')
    ax1.set_zlabel('Propability of Wavefunction')

    # Set from 0, 80 to 0, 1
    ax1.set_zlim3d(0, 1)
    if rotate==True:
        ax1.azim = k
    plt.savefig('visualization/figures/prob2D_%d.png' % (k))
    print('saving frame %d \n' % (k))


t2 = datetime.now()
delta = t2 - t1
print(delta)

#plt.show()
print('done plotting... \n')
print('running ffmpeg image to video shell script...')

#plt.show()

t2 = datetime.now()
delta = t2 - t1
print(delta)

rc = subprocess.call("ffmpeg -r 30 -f image2 -s 800x800 -i visualization/figures/prob2D_%01d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p output-prob2D.mp4",shell=True)
