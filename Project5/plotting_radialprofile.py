"""
RADIAL DENSITY PROFILE
A program that:
    
    1) Draws a histogram for the particle density as a function of radius,
    2) Plots the radial profile density,
    3) Compares the radial profile density with two radial profile density functions.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
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
    

def radial_profile(total_stars,time_step,final):
    # Input: Total stars, time step, whether to print final (True) or initial (False )distribution
    # Output: Plots and statistics.
    
    # Read file
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    t,index,m,x,y,z,vx,vy,vz = read_file(filename)
    
    if final == True:
        title = 'Final configuration'
        x = x[-total_stars:0]
        y = y[-total_stars:0]
        z = z[-total_stars:0]
    else:
        title = 'Intitial configuration'
        x = x[0:total_stars]
        y = y[0:total_stars]
        z = z[0:total_stars]
    
    # Calculate radial distances
    center = [0,0,0]
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
    
    # Print statistics
    print "Average radius = \t", np.average(r)
    print "Standard deviation = \t", np.std(r)
    
    # Plot histogram
    plt.figure()
    plt.hist(r,11)
    plt.xlabel('Radial distance [ly]')
    plt.ylabel('Frequency')
    plt.title(title,size=12)
    plt.show()
    
    """
    # Plot radial density profile
    r = r.astype(np.int)
    counts = np.bincount(r)
    radii = np.linspace(0,0.5,len(counts))
    
    plt.figure()
    plt.plot(radii,counts,label='Numerical')
    plt.xlabel('Radial distance [ly]',size=12)
    plt.ylabel('Radial density',size=12)
    plt.title('Radial density profile',size=12)
    
    # Make simple fit
    n_0 = 20
    r_0 = 1
    density_fit_simple = simple_fit(n_0,r_0,radii)
    
    # Make NFW-it
    rho_0 = 20
    r_0 = 3
    density_fit_NFW = NFW_profile(rho_0,r_0,radii)
    
    # Plot fits
    plt.plot(radii,density_fit_simple,label='Simple')
    plt.plot(radii,density_fit_NFW,label='NFW')
    plt.xlabel('Radius [ly]',size=12)
    plt.ylabel('Number density',size=12)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    """
    
    return


def simple_fit(n_0,r_0,r):
    return n_0/(1 + (r/r_0)**4)


def NFW_profile(rho_0,r_0,r):
    return rho_0/(r/r_0)/(1 + r/r_0)**2
    

def main(argv):
    
    total_stars = 100
    time_step = 0.05
    integration_points = 100

    # Radial profile
    final = False # False = intitial distribution, True = final distribution
    radial_profile(total_stars,time_step,final)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    