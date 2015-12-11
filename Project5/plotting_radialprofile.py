"""
RADIAL DENSITY PROFILE
A program that:
    
    1) Draws a histogram showing the radial distribution,
    2) Plots the radial profile density,
    3) Compares the radial profile density with a radial profile density function,
    4) Calculates average radial distance and standard deviation.
    
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
    

def radial_profile(total_stars,time_step,final,integration_points):
    # Input: Total stars, time step, whether to print final (True) or initial (False) distribution
    # Output: Plots and statistics.
    
    # Read file
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (total_stars,time_step)
    t,index,m,x,y,z,vx,vy,vz = read_file(filename)
    
    if final == True:
        title = 'Final configuration, N = %d, time step = %.3f' % (total_stars,time_step)
        x = x[total_stars*integration_points:-1]
        y = y[total_stars*integration_points:-1]
        z = z[total_stars*integration_points:-1]
        t = t[total_stars*integration_points:-1]
    else:
        title = 'Intitial configuration, N = %d, time step = %.3f' % (total_stars,time_step)
        x = x[0:total_stars]
        y = y[0:total_stars]
        z = z[0:total_stars]
        t = t[0:total_stars]
    
    # Calculate radial distances
    center = [0,0,0]
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        
    # Print radial statistics
    mu = np.average(r)
    stdev = np.std(r)
    
    print "Average radius = \t", mu
    print "Standard deviation = \t", stdev
    
    
    radial_distance = np.zeros(10)
    
    # Calculate the number of objects within given radii, in steps of 2 ly
    for i in range(total_stars-1):
        if r[i] <= 2:
            radial_distance[0] += 1
            
        elif r[i] > 2 and r[i] <= 4:
            radial_distance[1] += 1
            
        elif r[i] > 4 and r[i] <= 6:
            radial_distance[2] += 1
            
        elif r[i] > 6 and r[i] <= 8:
            radial_distance[3] += 1
            
        elif r[i] > 8 and r[i] <= 10:
            radial_distance[4] += 1
            
        elif r[i] > 10 and r[i] <= 12:
            radial_distance[5] += 1
        
        elif r[i] > 12 and r[i] <= 14:
            radial_distance[6] += 1
        
        elif r[i] > 14 and r[i] <= 16:
            radial_distance[7] += 1
        
        elif r[i] > 16 and r[i] <= 18:
            radial_distance[8] += 1
        
        elif r[i] > 18 and r[i] <=20:
            radial_distance[9] += 1

    
    # Calculate volume of shells, thickness 2 ly each
    shells = np.zeros(10)
    
    for i in range(10):
        shells[i] = volume((i+1)*2) - volume(i*2)
    
    # Calculate number densities
    no_density = radial_distance/shells
    
    # Calculate simple fit
    radii = np.linspace(0.001,20,len(no_density))
    radii_hires = np.linspace(0.001,20,1000)
    
    n_0 = no_density[0]
    r_0 = np.asarray([2,4,5,6])
    
    # Plot number density 
    plt.figure()
    plt.plot(radii/total_stars**(-1./3),no_density/(total_stars**2),'o',label='Simulation')
    plt.plot(radii_hires/total_stars**(-1./3),simple_fit(n_0,r_0[0],radii_hires)/(total_stars**2),label=r'Simple fit, $r_0 = 2$')
    plt.plot(radii_hires/total_stars**(-1./3),simple_fit(n_0,r_0[1],radii_hires)/(total_stars**2),label=r'Simple fit, $r_0 = 4$')
    plt.plot(radii_hires/total_stars**(-1./3),simple_fit(n_0,r_0[2],radii_hires)/(total_stars**2),label=r'Simple fit, $r_0 = 5$')
    plt.xlabel(r'$r/N^{-1/3}$')
    plt.ylabel(r'$n(r)/N^2$')
    plt.title(r'N = %d, time step = %.3f' % (total_stars,time_step),size=12)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=3,prop={'size':12})
    plt.show()
    
    # Plot histogram
    plt.figure()
    plt.hist(r,bins=20)
    plt.xlabel('Radial distance [ly]')
    plt.ylabel(r'Radial density $[N/$ly$^3]$')
    plt.xlim(0,20)
    plt.title(title,size=12)
    plt.show()
    
    # Plot average distance as function of N
    #N_values = [100,200,300,400]
    #averages = []
    
    #plt.figure()
    #plt.plot(N_values,averages)
    #plt.xlabel(r'$N$')
    #plt.ylabel(r'$\langle r \rangle$')
    #plt.show()
    
    return


def simple_fit(n_0,r_0,r):
    return n_0/(1 + (r/r_0)**4)


def NFW_profile(rho_0,r_0,r):
    return rho_0/(r/r_0)/(1 + r/r_0)**2

def volume(R):
    # Volume of sphere of radius R
    return 4*np.pi*R**3/3.0
    

def main(argv):
    
    total_stars = 100
    time_step = 0.001
    integration_points = 10000

    # Radial profile
    final = True # False = intitial distribution, True = final distribution
    radial_profile(total_stars,time_step,final,integration_points)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 
    