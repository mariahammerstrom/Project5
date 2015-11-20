"""
SCATTER PLOT
A program that creates a scatter plot in 3D of a cluster of stars, with stars
scaled according to their mass.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})


def read_file(filename):
    # Input: filename for file with structure [t,x,v_x,y,v_y,z,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    index = data[0]        # index
    m = data[1]            # mass
    
    # Positions
    x = data[2]            # position, x-direction
    y = data[3]            # position, y-direction
    z = data[4]            # position, z-direction
    
    # Velocities
    v_x = data[5]          # velocity, x-direction
    v_y = data[6]          # velocity, y-direction
    v_z = data[7]          # velocity, z-direction
    
    return index,m,x,y,z,v_x,v_y,v_z
    

def plot_scatter(N,x,y,z,mass):
    # Input: N = no. of stars, x,y,z = positions for each star, mass = mass of each star
    # Output: Scatter plot.
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Scale stars according to their mass
    for i in range(len(mass)):
        if mass[i] < 9:
            s = 2
        elif mass[i] > 9 and mass[i] < 11:
            s = 15
        elif mass[i] > 11:
            s = 30
        
        ax.scatter(x[i],y[i],z[i],s=s,c='k')

    ax.set_xlabel(r'$x$',size=14)
    ax.set_ylabel(r'$y$',size=14)
    ax.set_zlabel(r'$z$',size=14)
    
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    
    ax.set_title('Initial configuration',size=12)

    plt.show()



N = 100 # Number of stars

# Initial configuration
filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/galaxy_%d_20.0.txt' % N
index,m,x,y,z,v_x,v_y,v_z = read_file(filename)
    
plot_scatter(N,x,y,z,m)

# Final configuration

