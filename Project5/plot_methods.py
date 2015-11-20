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


def read_file(filename):
    # Input: filename for file with structure [t,x,y,z,v_x,v_y,v_z]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    
    # Positions
    x = data[1]            # position,x-direction
    y = data[2]            # position,y-direction
    z = data[3]            # position,z-direction
    
    # Velocities
    v_x = data[4]          # velocity, x-direction
    v_y = data[5]          # velocity, y-direction
    v_z = data[6]          # velocity, z-direction
    
    return t,x,y,z,v_x,v_y,v_z


def plot_time(N,time_step):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/VV_%.3f.txt' % time_step
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/RK4_%.3f.txt' % time_step
    
    t_verlet,x_verlet,y_verlet,z_verlet,v_x_verlet,v_y_verlet,v_z_verlet = read_file(filename_verlet)
    t_RK4,x_RK4,y_RK4,z_RK4,v_x_RK4,v_y_RK4,v_z_RK4 = read_file(filename_RK4)
    
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
    
    plt.figure(2)
    plt.title('Velocity, N = %d, time step = %.2f' % (N,time_step),size=12)
    plt.plot(t_verlet,v_x_verlet,label='Verlet')
    plt.plot(t_RK4,v_x_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$v_x$',size=14)
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


def main(argv):
    # Plot results as a function of time
    plot_time(N=2,time_step=0.1)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 