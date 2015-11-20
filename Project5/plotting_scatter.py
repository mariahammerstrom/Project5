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
    
    t = data[0]            # time
    x = data[1]            # position, x-direction
    v_x = data[2]          # velocity, x-direction
    y = data[3]            # position, y-direction
    v_y = data[4]          # velocity, y-direction
    z = data[5]            # position, z-direction
    v_z = data[6]          # velocity, z-direction
    m = data[7]            # mass
    
    return t,x,v_x,y,v_y,z,v_z,m
    

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
    
    ax.set_title('Initial configuration',size=12)

    plt.show()



N = 5 # Number of stars

# Initial configuration
x_coord = np.zeros(N)
y_coord = np.zeros(N)
z_coord = np.zeros(N)
#mass = np.zeros(N)
mass = np.array([8,10,12,15,5])

for i in range(N):
    #filename = 'stars_initial.txt'
    #t,x,v_x,y,v_y,z,v_z,m = read_file(filename)
    x_coord[i] = i
    y_coord[i] = i
    z_coord[i] = i
    #mass[i] = m

plot_scatter(N,x_coord,y_coord,z_coord,mass)

# Final configuration

