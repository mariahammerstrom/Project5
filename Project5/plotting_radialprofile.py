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

import matplotlib.mlab as mlab


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
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (total_stars,time_step)
    t,index,m,x,y,z,vx,vy,vz = read_file(filename)
    
    if final == True:
        title = 'Final configuration'
        x = x[-total_stars:-1]
        y = y[-total_stars:-1]
        z = z[-total_stars:-1]
        t = t[-total_stars:-1]
    else:
        title = 'Intitial configuration'
        x = x[0:total_stars]
        y = y[0:total_stars]
        z = z[0:total_stars]
        t = t[0:total_stars]
    
    # Calculate radial distances
    center = [0,0,0]
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
    
    volume = 4*np.pi*r**3/3.0;
    
    # Calculate radial densities
    r_rounded = np.around(r, decimals=0)
    counts = np.bincount(r_rounded.astype(int),weights=None,minlength=None)
    print r
    print r_rounded.astype(int)
    print counts
    
    # Print statistics
    mu = np.average(r)
    stdev = np.std(r)
    sigma = stdev**2
    
    print "Average radius = \t", mu
    print "Standard deviation = \t", stdev
    
    # Plot histogram
    bins = 60
    
    """
    plt.figure()
    plt.hist(r,bins=bins,normed=True)
    plt.xlabel('Radial distance [ly]')
    plt.ylabel('Frequency')
    plt.title(title,size=12)
    
    plt.figure()
    weights = 1./volume
    plt.hist(r,bins=bins,weights=weights)
    plt.xlabel('Radial distance [ly]')
    plt.ylabel(r'Radial density $[N/$ly$^3]$')
    plt.title(title,size=12)
    plt.show()
    
    plt.figure()
    plt.plot(r,'.')
    plt.title('r')
    plt.show()
    """
     

    # Make simple fit
    n_0 = total_stars**2
    r_0 = total_stars**(-1./3)
    radii = np.linspace(0,20,1000)
    density_fit_simple = simple_fit(n_0,r_0,radii)
    density_fit_NRW = NFW_profile(n_0,r_0,radii)

    # Plot fits
    plt.figure()
    plt.plot(radii/total_stars**(-1./3),density_fit_simple/(total_stars**2),label='Simple')
    plt.plot(radii/total_stars**(-1./3),density_fit_NRW/(total_stars**2),label='NRW')
    plt.xlabel('Radius [ly]',size=12)
    plt.ylabel('Number density',size=12)
    plt.legend(loc=1,prop={'size':12})
    plt.xlim(0,20)
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

    # the histogram of the data
    plt.figure()    
    weights = 1./volume/total_stars**(-1./3) 
    n, bins, patches = plt.hist(r, 60, normed=1, facecolor='green', alpha=0.75, weights=weights)
    
    # Plot
    plt.xlabel('Radial distance [ly]')
    plt.ylabel(r'Radial density $[N/$ly$^3]$')
    plt.xlim(0,20)
    plt.grid(True)

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
    plt.show()
    """
    
    return


def simple_fit(n_0,r_0,r):
    return n_0/(1 + (r/r_0)**4)


def NFW_profile(rho_0,r_0,r):
    return rho_0/(r/r_0)/(1 + r/r_0)**2
    

def main(argv):
    
    total_stars = 500
    time_step = 0.005
    integration_points = 100

    # Radial profile
    final = False # False = intitial distribution, True = final distribution
    radial_profile(total_stars,time_step,final)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    