"""
PLOT: N-BODY PROBLEM (3D)
A program that plots the position x and velocity v as a function of time t for N masses in a gravitational field
solved with 4th order Runge-Kutta (RK4) and Velocity-Verlet (VV), compared with the analytical solution.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


def read_file(filename,total_stars):
    # Input: filename for file with structure [t,star no.,m,x,y,z,vx,vy,vz]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    # Time step
    t = data[0]
    time_steps = t[::total_stars]
    
    # Masses
    mass = data[1]
    m = np.array(len(total_stars))
    
    for i in range(total_stars):
        m[i] = mass[i]
    
    # Positions
    x = data[2]            # position,x-direction
    y = data[3]            # position,y-direction
    z = data[4]            # position,z-direction
    
    # Velocities
    vx = data[5]
    vy = data[6]          # velocity, y-direction
    vz = data[7]          # velocity, z-direction
    
    return time_steps,m,x,y,z,vx,vy,vz

    

def plot_time(stars,N,time_step):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_VV_%d_%.2f.txt' % (N,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_RK4_%d_%.2f.txt' % (N,time_step)
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (stars,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.2f.txt' % (stars,time_step)
    
    t_verlet,x_verlet,y_verlet,z_verlet,v_x_verlet,v_y_verlet,v_z_verlet,r_verlet = read_file2(filename_verlet)
    t_RK4,x_RK4,y_RK4,z_RK4,v_x_RK4,v_y_RK4,v_z_RK4,r_RK4 = read_file2(filename_RK4)
    
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
    
    # plt.figure(4)
    # plt.title('Force, N = %d, time step = %.2f' % (N,time_step),size=12)
    # plt.plot(t_verlet,F_verlet,label='Verlet')
    # plt.plot(t_RK4,F_RK4,label='RK4')
    # plt.xlabel(r'$t$',size=14)
    # plt.ylabel(r'$F_x$',size=14)
    # plt.legend(loc=2,prop={'size':12})
    # plt.show()
    
    plt.figure(3)
    plt.title('Total energy, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,np.array(E_verlet)-0.5,label='Verlet')
    plt.plot(t_RK4,np.array(E_RK4)-0.5,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$E$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    return


def plot_orbits(stars,N,time_step):
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.2f.txt' % (stars,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_%d_%.2f.txt' % (stars,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.2f.txt' % (stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.2f.txt' % (stars,time_step)
    
    # t1_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet,F_verlet = read_file(filename_verlet)
    # t1_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4,F_RK4 = read_file(filename_RK4)
    
    t1_verlet,x1_verlet,y1_verlet,z1_verlet,v_x1_verlet,v_y1_verlet,v_z1_verlet,r1_verlet = read_file2(filename_verlet)
    t1_RK4,x1_RK4,y1_RK4,z1_RK4,v_x1_RK4,v_y1_RK4,v_z1_RK4,r1_RK4 = read_file2(filename_RK4)
    
    mpl.rcParams['legend.fontsize'] = 10
 
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    limit = -1 #70 for large time step
    ax.plot(x1_verlet[0:limit],y1_verlet[0:limit],z1_verlet[0:limit],'r', label='#1, VV')
    ax.plot(x1_RK4[0:limit],y1_RK4[0:limit],z1_RK4[0:limit], 'b',label='#1, RK4')
    
    #ax.plot(x2_verlet,y2_verlet,z2_verlet,'g', label='#2, VV')
    #ax.plot(x2_RK4,y2_RK4,z2_RK4, 'y',label='#2, RK4')
    ax.legend()

    plt.show()

    return
    
def main(argv):
    # Plot results as a function of time
    #plot_time(stars=2,N=10000,time_step=.1)
    #plot_orbits(stars=2,N=1000,time_step=1)
    
    total_stars = 2
    time_step = 1
    filename = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/a_cluster_VV_%d_%.2f.txt' % (total_stars,time_step)
    time_steps,m,x,y,z,vx,vy,vz = read_file(filename,total_stars)
    
    print time_steps
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 