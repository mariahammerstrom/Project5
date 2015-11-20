"""
SCATTER PLOT
A program that creates a scatter plot in 3D of a cluster of stars, with stars
scaled according to their mass.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
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
    

def plot_scatter(N,x,y,z,mass,title):
    # Input: N = no. of stars, x,y,z = positions for each star, mass = mass of each star
    # Output: Scatter plot.
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Scale stars according to their mass
    for i in range(len(mass)):
        if mass[i] < 9:
            #s = 2
            c = 'c'
        elif mass[i] > 9 and mass[i] < 11:
            #s = 15
            c = 'y'
        elif mass[i] > 11:
            #s = 30
            c = 'm'
        
        ax.scatter(x[i],y[i],z[i],s=30,c=c)

    ax.set_xlabel(r'$x$',size=14)
    ax.set_ylabel(r'$y$',size=14)
    ax.set_zlabel(r'$z$',size=14)
    
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    
    ax.set_title(title,size=12)

    plt.show()



def main(argv):
    N = 100 # Number of stars

    # Read data
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/galaxy_%d_20.0.txt' % N
    index,m,x,y,z,v_x,v_y,v_z = read_file(filename)
    
    # Initial configuration
    plot_scatter(N,x[0:100],y[0:100],z[0:100],m[0:100],'Initial condition, t = 0')

    # Final configuration
    #print index[-100:]
    #plot_scatter(N,x[-100:],y[-100:],z[-100:],m[-100:],'Final condition, t = %.2f' % index[-1])
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    