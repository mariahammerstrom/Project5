"""
PLOT: 2-BODY PROBLEM (3D)
A program that plots the position x and velocity v as a function of time t for two masses in a gravitational field
solved with 4th order Runge-Kutta (RK4) and Velocity-Verlet (VV), compared with the analytical solution.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from math import pi

G = 4*pi**2


def read_file(filename):
    # Input: filename for file with structure [t,x,y,z,v_x,v_y,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    
    # Positions
    x = data[3]            # position,x-direction
    v_x = data[4]

    y = data[5]            # position,y-direction
    v_y = data[6]          # velocity, y-direction
    
    z = data[7]            # position,z-direction
    v_z = data[8]          # velocity, z-direction
    
    # Force
    F = data[9]
    
    r = np.sqrt(x**2 + y**2 + z**2)
    
    return t,x,y,z,v_x,v_y,v_z,r,F

def read_file2(filename):
    # Input: filename for file with structure [t,x,y,z,v_x,v_y,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    m = data[2]			   # mass
    
    # Positions
    x = data[3]            # position,x-direction
    y = data[4]            # position,y-direction
    z = data[5]            # position,z-direction
    
    # Velocities
    v_x = data[6]          # velocity, x-direction
    v_y = data[7]          # velocity, y-direction
    v_z = data[8]          # velocity, z-direction
    
    r = np.sqrt(x**2 + y**2 + z**2)
    
    return t,m,x,y,z,v_x,v_y,v_z,r
    

def plot_time(stars,N,time_step):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (stars,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_%d_%.3f.txt' % (stars,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.3f.txt' % (stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.3f.txt' % (stars,time_step)
    
    t_verlet,m_verlet,x_verlet,y_verlet,z_verlet,v_x_verlet,v_y_verlet,v_z_verlet,r_verlet = read_file2(filename_verlet)
    t_RK4,m_RK4,x_RK4,y_RK4,z_RK4,v_x_RK4,v_y_RK4,v_z_RK4,r_RK4 = read_file2(filename_RK4)
    
    # Calculate energies
    E_verlet = 0.5*np.array(v_x_verlet)**2 - G*np.array(m_verlet)**2/np.array(r_verlet)
    E_RK4 = 0.5*np.array(v_x_RK4)**2 - G*np.array(m_RK4)**2/np.array(r_verlet)
    
    # Make plots
    plt.figure()
    plt.title('Position, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,x_verlet,'r',label='VV')
    plt.plot(t_RK4,x_RK4,'b',label='RK4')
    plt.xlabel(r'$t$ $\mathrm{[year]}$',size=14)
    plt.ylabel(r'$x$ $\mathrm{[AU]}$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    """
    plt.figure()
    plt.title('Radial position, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,r_verlet,label='Verlet')
    plt.plot(t_RK4,r_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$r$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    """
    
    plt.figure()
    plt.title('Velocity, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,v_x_verlet,'r',label='VV')
    plt.plot(t_RK4,v_x_RK4,'b',label='RK4')
    plt.xlabel(r'$t$ $\mathrm{[year]}$',size=14)
    plt.ylabel(r'$v_x$ $\mathrm{[AU/year]}$',size=14)
    plt.legend(loc=2,prop={'size':12})
    plt.show()
    
    # plt.figure()
    # plt.title('Force, N = %d, time step = %.2f' % (N,time_step),size=12)
    # plt.plot(t_verlet,F_verlet,label='Verlet')
    # plt.plot(t_RK4,F_RK4,label='RK4')
    # plt.xlabel(r'$t$',size=14)
    # plt.ylabel(r'$F_x$',size=14)
    # plt.legend(loc=2,prop={'size':12})
    # plt.show()
    
    plt.figure()
    plt.title('Total energy, N = %d, time step = %.3f' % (N,time_step),size=12)
    plt.plot(t_verlet,np.array(E_verlet),'r',label='VV')
    plt.plot(t_RK4,np.array(E_RK4),'b',label='RK4')
    plt.xlabel(r'$t$ $\mathrm{[year]}$',size=14)
    plt.ylabel(r'$E$ $[M_{\odot}\mathrm{AU}^2\mathrm{/years}^2]$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    return


def plot_orbits(stars,N,time_step):
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (stars,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_%d_%.3f.txt' % (stars,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.2f.txt' % (stars,time_step)
    
    # t1_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet,F_verlet = read_file(filename_verlet)
    # t1_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4,F_RK4 = read_file(filename_RK4)
    
    t1_verlet,m_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet = read_file2(filename_verlet)
    t1_RK4,m_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4 = read_file2(filename_RK4)
    
    mpl.rcParams['legend.fontsize'] = 10
 
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    limit = -1 #70 for large time step
    ax.plot(x1_verlet[0:limit],y1_verlet[0:limit],z1_verlet[0:limit],'r', label='VV')
    ax.plot(x1_RK4[0:limit],y1_RK4[0:limit],z1_RK4[0:limit], 'b',label='RK4')
    
    ax.set_xlabel(r'$x$ $\mathrm{[AU]}$')
    ax.set_ylabel(r'$y$ $\mathrm{[AU]}$')
    ax.set_zlabel(r'$z$ $\mathrm{[AU]}$')
    ax.legend()

    plt.show()

    return


def plot_orbits2D(stars,N,time_step):
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (stars,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_%d_%.3f.txt' % (stars,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.2f.txt' % (stars,time_step)
    
    # t1_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet,F_verlet = read_file(filename_verlet)
    # t1_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4,F_RK4 = read_file(filename_RK4)
    
    t1_verlet,m_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet = read_file2(filename_verlet)
    t1_RK4,m_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4 = read_file2(filename_RK4)
    
    mpl.rcParams['legend.fontsize'] = 10
 
    plt.figure()
    limit = -1 #70 for large time step
    plt.plot(x1_verlet[0:limit],y1_verlet[0:limit],'r', label='VV')
    plt.plot(x1_RK4[0:limit],y1_RK4[0:limit], 'b',label='RK4')
    plt.title('Earth, time step = %.3f' % time_step,size=12)
    
    plt.xlabel(r'$x$ $\mathrm{[AU]}$')
    plt.ylabel(r'$y$ $\mathrm{[AU]}$')
    plt.xlim(-1.2,1.2)
    plt.ylim(-1.2,1.2)
    plt.legend()

    plt.show()

    return
    
def main(argv):
    # Plot results as a function of time
    #plot_time(stars=2,N=10,time_step=0.5)
    #plot_orbits(stars=2,N=10000,time_step=0.005)
    plot_orbits2D(stars=2,N=10000,time_step=0.005)
    
    print 6.3*(1.49e8)/(365*24*60*60) # [km/s], should be approx. 29.78 km/s
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 