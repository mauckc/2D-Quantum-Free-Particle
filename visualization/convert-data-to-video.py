# Visualizing a 2D scalar field
# matplotlib and NumPy offer some interesting mechanisms that make the visualization of a 2D scalar field convenient.
# In this recipe, we show a very simple way to visualize a 2D scalar field.

# Setup Dependecies:
#
# pip install matplotlib
# sudo apt install python-tk

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import subprocess

# Set dimension of our 2D grid
n = 128
noutputs = 100
# Load data files from slices output directory
DataIn = np.loadtxt('../slices/slices_fields_0.dat')

# Create array for storing each field dat file data
AllSlicesData = np.zeros((noutputs, n*n, 5))
Allprob = np.zeros((noutputs, n*n))
# For number of outputs user variable load a successive slice
for x in xrange(noutputs):
    AllSlicesData[x] = np.loadtxt('../slices/slices_fields_%d.dat' % (x))
    Allprob[x] = AllSlicesData[x,::,4]

# Interpret All slices data fields
indexA = AllSlicesData[-1,::,0]
indexB = AllSlicesData[-1,::,1]
phi = AllSlicesData[-1,::,2]
chi = AllSlicesData[-1,::,3]
prob = AllSlicesData[-1,::,4]

for p in xrange(10):
    print Allprob[p]

# Interpret single slice data fields
#indexA = DataIn[::,0]
#indexB = DataIn[::,1]
#phi = DataIn[::,2]
#chi = DataIn[::,3]
#prob = DataIn[::,4]

# Create figure for plot
fig = plt.figure(figsize=(8,8))

# Populate 1-Ds arrays with length n of floating point values between 0 and 1 for our bounds
x = np.linspace(0., 1., n)
y = np.linspace(0., 1., n)
#The numpy.meshgrid() function generates the samples from an explicit 2D function.
X, Y = np.meshgrid(x, y)
Z = np.zeros((n, n))

for x in xrange(noutputs):
    # Must interpret our indexed array and save it to the new data type
    for j in xrange(0, n-1):
        for i in xrange(0, n-1):
            Z[i][j] = Allprob[x,i+j*(n-1)]
    #Then, pyplot.pcolormesh() is used to display the function, as shown in the following code:
    plt.pcolormesh(X, Y, Z, cmap = cm.gray)
    plt.savefig('figures/prob2D_%d.png' % (x))
    print('saving frame %d \n' % (x))

plt.show()

rc = subprocess.call("./compilepngs2video.sh",shell=True)

# To compile the images into a video use ffmpeg
# ffmpeg can be installed on ubuntu using "sudo apt-get install ffmpeg"
# ffmpeg -r 30 -f image2 -s 800x800 -i figures/prob2D_%01d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4

# This wasn't working because my newly installed ubuntu 16.04 did not have the avi or mp4 codecs linked to ffmpeg
# create video out of the images
#import cv2

#out = cv2.VideoWriter('output.avi', cv2.VideoWriter_fourcc('M','J','P','G'), 20.0, (800,800))
# create array to store all images *might need to add number for RGBA if cmap wasnt grey
#image = np.zeros((noutputs,800,800,3))
#for x in xrange (noutputs):
    #image[x] = cv2.imread('figures/prob2D_%d.png' % (x) )

#for x in xrange (noutputs):
#    out.write(image[x])

#out.release()
