"""
SCATTER PLOT
A program that creates scatter plot(s) in 3D of a cluster of stars, either
at a given time (initial or final) or as a function of time for the whole calculation.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import mpl_toolkits.mplot3d
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})



def read_file(filename):
    # Input: filename for file with structure [t,index,m,x,y,z,vx,vy,vz]
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
    

def plot_scatter(total_stars,x,y,z,title,show):
    # Input: total_stars = total no. of stars, x,y,z = positions for each star, title = plot title,
    # show = True if plot should be shown, False if plot should not be shown.
    # Output: Scatter plot.
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    s = ax.scatter(x,y,z,s=30,color='purple',marker='o')
    s.set_edgecolors = s.set_facecolors = lambda *args:None

    ax.set_xlabel(r'$x$',size=14)
    ax.set_ylabel(r'$y$',size=14)
    ax.set_zlabel(r'$z$',size=14)
    
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    ax.set_zlim(-20, 20)
    
    ax.set_title(title,size=12)
    
    if show == True:
        plt.show()


def scatter_plot_series(total_stars,time_step,integration_points,x,y,z,t):
    # Input: total number of stars, time step, integration points, positions (x,y,z), time array
    # Output: scatter plots as a function of time for each integration point.
    
    # Delete previously stored image files
    delete_files(total_stars)
    
    # Plot scatter plot as time evolves for the system and save plots to file.
    for i in range(integration_points):
        plot_scatter(total_stars,x[i*total_stars:(i+1)*total_stars],y[i*total_stars:(i+1)*total_stars],z[i*total_stars:(i+1)*total_stars],'t = %.2f' % t[i*total_stars],False)
        plt.savefig('Movie/Image_%d_%.2f.png' % (total_stars,t[i*total_stars]))
        
    return


def delete_files(total_stars):
    # Function to delete previously stored files
    
    import glob,os
    for f in glob.glob('Movie/Image_%d_*.png' % total_stars):
        os.remove(f)
    return



def main(argv):
    total_stars = 100
    time_step = 0.05
    integration_points = 100
    
    # Read data
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    #filename = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    t,index,m,x,y,z,v_x,v_y,v_z = read_file(filename)
    
    # Initial configuration
    plot_scatter(total_stars,x[0:total_stars],y[0:total_stars],z[0:total_stars],'Initial condition, t = 0',True)

    # Final configuration
    #plot_scatter(total_stars,x[-total_stars:],y[-total_stars:],z[-total_stars:],'Final condition, t = %.2f' % t[-1],True)
    
    # Series of scatter plots as a function of time
    #scatter_plot_series(total_stars,time_step,integration_points,x,y,z,t)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    