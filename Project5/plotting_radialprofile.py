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
    

def radial_profile(x,y,z,center):
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
    
    # Print statistics
    print "Average radius = \t", np.average(r)
    print "Standard deviation = \t", np.std(r)
    
    # Plot histogram
    plt.figure()
    plt.hist(r,21)
    plt.xlabel('Radial distance')
    plt.ylabel('Frequency')
    plt.show()
    
    # Plot radial density profile
    r = r.astype(np.int)
    counts = np.bincount(r)
    radii = np.linspace(0,20,len(counts))
    
    plt.figure()
    plt.plot(radii,counts,label='Numerical')
    plt.xlabel('Radius [ly]',size=12)
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
    
    return


def simple_fit(n_0,r_0,r):
    return n_0/(1 + (r/r_0)**4)


def NFW_profile(rho_0,r_0,r):
    return rho_0/(r/r_0)/(1 + r/r_0)**2
    

def main(argv):
    
    # Read file
    N = 100 # Number of stars
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/galaxy_%d_20.0.txt' % N
    index,m,x,y,z,v_x,v_y,v_z = read_file(filename)
    
    # Radial profile
    center = [0,0,0]
    radial_profile(x,y,z,center)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    