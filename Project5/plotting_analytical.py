# PROJECT 5a)

"""
A program that plots the position x and velocity v as a function of time t for a box on a spring
solved with 4th order Runge-Kutta (RK4) and Velocity-Verlet (VV), compared with the analytical solution.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})



def read_file(filename):
    # Input: filename for file with structure [t,x,v]
    # Output: arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    t = data[0]            # time
    x = data[1]            # position
    v = data[2]            # velocity
    
    return t,x,v
    

def plot_time(time_step):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/Verlet_%.1f.txt' % time_step
    t_verlet,x_verlet,v_verlet = read_file(filename_verlet)
    
    E_verlet = 0.5*np.array(v_verlet)**2 + 0.5*np.array(x_verlet)**2
    
    # Make plots
    plt.figure(1)
    plt.title('Comparing methods (position), time step = %.1f' % time_step,size=12)
    plt.plot(t_verlet,np.cos(t_verlet),label='Analytic')
    plt.plot(t_verlet,x_verlet,label='Verlet')
    #plt.plot(t_RK4,x_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$x$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    plt.figure(2)
    plt.title('Comparing methods (velocity), time step = %.1f' % time_step,size=12)
    plt.plot(t_verlet,-np.sin(t_verlet),label='Analytic')
    plt.plot(t_verlet,v_verlet,label='Verlet')
    #plt.plot(t_RK4,v_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$v_x$',size=14)
    plt.legend(loc=2,prop={'size':12})
    plt.show()
    
    #plt.figure(3)
    #plt.title('Total energy',size=12)
    #plt.plot(t_verlet,E_verlet)
    #plt.xlabel(r'$t$',size=14)
    #plt.ylabel(r'$E$',size=14)
    #plt.show()
    
    return


def plot_timestep(time_step_list):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution for different time steps.
    
    # Define plots
    plt.figure(1)
    plt.title('Comparing methods (position)',size=12)
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$x$',size=14)
    
    plt.figure(2)
    plt.title('Comparing methods (velocity)',size=12)
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$v_x$',size=14)
    
    
    for i in range(len(time_step_list)):
        filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/Verlet_%.1f.txt' % time_step_list[i]
        t_verlet,x_verlet,v_verlet = read_file(filename_verlet)
        
        plt.figure(1)
        plt.plot(t_verlet,np.cos(t_verlet),label='Analytic, dt = %.1f' % time_step_list[i])
        plt.plot(t_verlet,x_verlet,label='Verlet, dt = %.1f' % time_step_list[i])
        #plt.plot(t_RK4,x_RK4,label='RK4')
        plt.legend(loc=1,prop={'size':12})
        
        plt.figure(2)
        plt.plot(t_verlet,-np.sin(t_verlet),label='Analytic, dt = %.1f' % time_step_list[i])
        plt.plot(t_verlet,v_verlet,label='Verlet, dt = %.1f' % time_step_list[i])
        #plt.plot(t_RK4,v_RK4,label='RK4')
        plt.legend(loc=2,prop={'size':12})
    
    plt.show()
    
    return 0


# Plot results as a function of time
time_step = 0.2
plot_time(time_step)

# Plot results for different time steps
#time_step_list = [0.2]
#plot_timestep(time_step_list)
