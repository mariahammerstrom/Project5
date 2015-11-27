"""
A program that plots the position x and velocity v as a function of time t for a box on a spring
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
    
    #t = data[0]            # time
    #x = data[1]            # position,x-direction
    #v_x = data[2]          # velocity, x-direction

    t = data[0]
    x = data[3]
    v_x = data[4]
    
    return t,x,v_x
    

def plot_time(time_step):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/VV_analytic_%.3f.txt' % time_step
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/RK4_analytic_%.3f.txt' % time_step
    filename_RK4 = 'C:/Users/mariefoss/Documents/GitHub/Project5/build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/analytic_RK4.txt'

    #t_verlet,x_verlet,v_verlet = read_file(filename_verlet)
    t_RK4,x_RK4,v_RK4 = read_file(filename_RK4)
    
    # Calculate energies
    #E_verlet = 0.5*np.array(v_verlet)**2 + 0.5*np.array(x_verlet)**2
    E_RK4 = 0.5*np.array(v_RK4)**2 + 0.5*np.array(x_RK4)**2
    E_analytic = 0.5*np.array(-np.sin(t_RK4))**2 + 0.5*np.array(np.cos(t_RK4))**2
    
    # Make plots
    plt.figure(1)
    plt.title('Comparing methods (position), time step = %.2f' % time_step,size=12)
    plt.plot(t_RK4,np.cos(t_RK4),label='Analytic')
    #plt.plot(t_verlet,x_verlet,label='Verlet')
    plt.plot(t_RK4,x_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$x$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    plt.figure(2)
    plt.title('Comparing methods (velocity), time step = %.2f' % time_step,size=12)
    plt.plot(t_RK4,-np.sin(t_RK4),label='Analytic')
    #plt.plot(t_verlet,v_verlet,label='Verlet')
    plt.plot(t_RK4,v_RK4,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$v_x$',size=14)
    plt.legend(loc=2,prop={'size':12})
    plt.show()
    
    plt.figure(3)
    plt.title('Total energy, time step = %.2f' % time_step,size=12)
    plt.plot(t_RK4,np.array(E_analytic)-0.5,label='Analytic')
    #plt.plot(t_verlet,np.array(E_verlet)-0.5,label='Verlet')
    plt.plot(t_RK4,np.array(E_RK4)-0.5,label='RK4')
    plt.xlabel(r'$t$',size=14)
    plt.ylabel(r'$E$',size=14)
    plt.legend(loc=1,prop={'size':12})
    plt.show()
    
    return


def main(argv):
    # Plot results as a function of time
    plot_time(time_step=0.1)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 