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


def read_file(filename):
    # Input: filename for file with structure [t,x,y,z,v_x,v_y,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    
    # Positions
    x = data[1]            # position,x-direction
    v_x = data[2]

    y = data[3]            # position,y-direction
    v_y = data[4]          # velocity, y-direction
    
    z = data[5]            # position,z-direction
    v_z = data[6]          # velocity, z-direction
    
    # Force
    F = data[7]
    
    r = np.sqrt(x**2 + y**2 + z**2)
    
    return t,x,y,z,v_x,v_y,v_z,r,F

def read_file2(filename):
    # Input: filename for file with structure [t,x,y,z,v_x,v_y,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    
    # Positions
    x = data[7]            # position,x-direction
    y = data[9]            # position,y-direction
    z = data[11]            # position,z-direction
    
    # Velocities
    v_x = data[8]          # velocity, x-direction
    v_y = data[10]          # velocity, y-direction
    v_z = data[12]          # velocity, z-direction
    
    r = np.sqrt(x**2 + y**2 + z**2)
    
    return t,x,y,z,v_x,v_y,v_z,r
    

def plot_time(N,time_step):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_VV_%d_%.2f.txt' % (N,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_RK4_%d_%.2f.txt' % (N,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.2f.txt' % (stars,time_step)
    
    t_verlet,x_verlet,y_verlet,z_verlet,v_x_verlet,v_y_verlet,v_z_verlet,r_verlet,F_verlet = read_file(filename_verlet)
    t_RK4,x_RK4,y_RK4,z_RK4,v_x_RK4,v_y_RK4,v_z_RK4,r_RK4,F_RK4 = read_file(filename_RK4)
    
    # Calculate energies
    E_verlet = 0.5*np.array(v_x_verlet)**2 + 0.5*np.array(x_verlet)**2
    E_RK4 = 0.5*np.array(v_x_RK4)**2 + 0.5*np.array(x_RK4)**2
    
    # Make plots
    plt.figure(1)
    plt.title('Position, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,x_verlet,label='Verlet')
    plt.plot(t_RK4,x_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$x$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    plt.figure(3)
    plt.title('Radial position, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,r_verlet,label='Verlet')
    plt.plot(t_RK4,r_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$r$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    plt.figure(2)
    plt.title('Velocity, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,v_x_verlet,label='Verlet')
    plt.plot(t_RK4,v_x_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$v_x$',size=14)
    plt.legend(loc=2,prop={'size':12})
    plt.show()
    
    plt.figure(4)
    plt.title('Force, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,F_verlet,label='Verlet')
    plt.plot(t_RK4,F_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$F_x$',size=14)
    plt.legend(loc=2,prop={'size':12})
    plt.show()
    
    plt.figure(3)
    plt.title('Total energy, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,np.array(E_verlet)-0.5,label='Verlet')
    plt.plot(t_RK4,np.array(E_RK4)-0.5,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$E$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    return
    
    
def plot_timestep(N,time_step_list):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution for different time steps.
    
    # Define plots
    plt.figure(1)
    plt.title('Position, N = $d' % N,size=12)
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$x$',size=14)
    
    plt.figure(2)
    plt.title('Velocity, N = $d' % N,size=12)
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$v_x$',size=14)
    
    
    for i in range(len(time_step_list)):
        filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/Verlet_%.1f.txt' % time_step_list[i]
        t_verlet,x_verlet,y_verlet,z_verlet,v_x_verlet,v_y_verlet,v_z_verlet = read_file(filename_verlet)
        
        plt.figure(1)
        plt.plot(t_verlet,x_verlet,label='Verlet, dt = %.1f' % time_step_list[i])
        #plt.plot(t_RK4,x_RK4,label='RK4')
        plt.legend(loc=1,prop={'size':12})
        
        plt.figure(2)
        plt.plot(t_verlet,v_x_verlet,label='Verlet, dt = %.1f' % time_step_list[i])
        #plt.plot(t_RK4,v_x_RK4,label='RK4')
        plt.legend(loc=2,prop={'size':12})
    
    plt.show()
    
    return 0


def plot_orbits(N,time_step):
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_VV_%d_%.2f.txt' % (N,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_RK4_%d_%.2f.txt' % (N,time_step)
    
    t1_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet,F_verlet = read_file(filename_verlet)
    t1_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4,F_RK4 = read_file(filename_RK4)
    
    #t2_verlet,x2_verlet,y2_verlet,z2_verlet,v_x2_verlet,v_y2_verlet,v_z2_verlet,r2_verlet = read_file2(filename_verlet)
    #t2_RK4,x2_RK4,y2_RK4,z2_RK4,v_x2_RK4,v_y2_RK4,v_z2_RK4,r2_RK4 = read_file2(filename_RK4)
    
    mpl.rcParams['legend.fontsize'] = 10
 
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x1_verlet,y1_verlet,z1_verlet,'r', label='#1, VV')
    ax.plot(x1_RK4,y1_RK4,z1_RK4, 'b',label='#1, RK4')
    
    #ax.plot(x2_verlet,y2_verlet,z2_verlet,'g', label='#2, VV')
    #ax.plot(x2_RK4,y2_RK4,z2_RK4, 'y',label='#2, RK4')
    ax.legend()

    plt.show()

    return
    
def main(argv):
    # Plot results as a function of time
    plot_time(N=10000,time_step=1.0)
    #plot_orbits(N=10000,time_step=1.0)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 