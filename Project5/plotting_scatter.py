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
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D


def read_file(filename):
    # Input: filename for file with structure [t,index,x,v_x,y,v_y,z,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    index = data[1]        # index
    m = data[2]            # mass
    
    # Positions
    x = data[3]            # position, x-direction
    y = data[4]            # position, y-direction
    z = data[5]            # position, z-direction
    
    # Velocities
    v_x = data[6]          # velocity, x-direction
    v_y = data[7]          # velocity, y-direction
    v_z = data[8]          # velocity, z-direction
    
    return t,index,m,x,y,z,v_x,v_y,v_z
    

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
    ax.set_zlim(-20, 20)
    
    ax.set_title(title,size=12)

    plt.show()



def main(argv):
    total_stars = 100
    time_step = 0.003
    integration_points = 1000
    
    #N = 100 # Number of stars

    # Read data
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    filename = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    t,index,m,x,y,z,v_x,v_y,v_z = read_file(filename)
    
    # Initial configuration
    plot_scatter(total_stars,x[0:total_stars],y[0:total_stars],z[0:total_stars],m[0:total_stars],'Initial condition, t = 0')

    # Final configuration
    #print index[-100:]
    plot_scatter(total_stars,x[-100:],y[-100:],z[-100:],m[-100:],'Final condition, t = %.2f' % t[-1])
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    