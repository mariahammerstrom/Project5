"""
SCATTER PLOT MOVIE!
A program that creates a scatter plot in 3D of a cluster of stars, with stars
scaled according to their mass.
"""

import sys,glob,os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D

from scitools.all import *

def read_file2(filename,star_number,total_stars):
    # Input: filename for file with structure [t,star no.,m,x,y,z,vx,vy,vz]
    # Output: arrays of data values for a given star
    
    # Read data
    data = np.loadtxt(filename,unpack=True) 
    
    # Time steps
    t = data[0]
    time = t[::total_stars]

    # Mass
    mass = data[2]
    m = mass[star_number]
    
    # Positions
    x = data[3]            # position,x-direction
    y = data[4]            # position,y-direction
    z = data[5]            # position,z-direction

    # Velocities
    vx = data[6]          # velocity, x-direction
    vy = data[7]          # velocity, y-direction
    vz = data[8]          # velocity, z-direction
    
    # Extract data for given star number
    list_range = len(time)
    x_star = np.zeros(list_range)
    y_star = np.zeros(list_range)
    z_star = np.zeros(list_range)
    vx_star = np.zeros(list_range)
    vy_star = np.zeros(list_range)
    vz_star = np.zeros(list_range)
    
    for i in range(len(time)):
            x_star[i] = x[i*total_stars + star_number]
            y_star[i] = y[i*total_stars + star_number]
            z_star[i] = z[i*total_stars + star_number]
            vx_star[i] = vx[i*total_stars + star_number]
            vy_star[i] = vy[i*total_stars + star_number]
            vz_star[i] = vz[i*total_stars + star_number]
                
    return time,m,x_star,y_star,z_star,vx_star,vy_star,vz_star
 
 
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


def plot_scatter(total_stars,x,y,z,title):
    # Input: N = no. of stars, x,y,z = positions for each star, mass = mass of each star
    # Output: Scatter plot.
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x,y,z,s=30,color='purple',marker='o')
    
    ax.set_xlabel(r'$x$',size=14)
    ax.set_ylabel(r'$y$',size=14)
    ax.set_zlabel(r'$z$',size=14)
    
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    ax.set_zlim(-20, 20)
    
    ax.set_title(title,size=12)
    
    return


def scatter_plot_series(total_stars,time_step,integration_points):
    
    for f in glob.glob('Movie/Image_*.eps'):
        os.remove(f)
    
    # Get data
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    #filename = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    t,index,m,x,y,z,v_x,v_y,v_z = read_file(filename)
    
    for i in range(integration_points):
        plot_scatter(total_stars,x[i*total_stars:(i+1)*total_stars],y[i*total_stars:(i+1)*total_stars],z[i*total_stars:(i+1)*total_stars],'t = %.2f' % t[i*total_stars])
        plt.savefig('Movie/Image_%d_%.2f.png' % (total_stars,t[i*total_stars]))
        
    return
    

def make_movie(total_stars,time_step):
    moviename = '/Movie/movie_%d.gif' % (total_stars)
    movie('Movie/Image_*.png',output_file=moviename)
    return 0


def main(argv):
    total_stars = 5
    time_step = 1.0
    integration_points = 9
    
    #scatter_plot_series(total_stars,time_step,integration_points)
    make_movie(total_stars,time_step)
    
    

    	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    