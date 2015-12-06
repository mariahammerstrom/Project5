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


def read_file(filename,star_number,total_stars):
    # Input: filename for file with structure [t,star no.,m,x,y,z,vx,vy,vz]
    # Output: arrays of data values for a given star
    
    # Read data
    data = np.loadtxt(filename,unpack=True) 
    
    # Time steps
    t = data[0]
    time = t[::total_stars]

    # Mass
    mass = data[2]
    m = mass[star_number]
    
    # Positions
    x = data[3]            # position,x-direction
    y = data[4]            # position,y-direction
    z = data[5]            # position,z-direction

    # Velocities
    vx = data[6]          # velocity, x-direction
    vy = data[7]          # velocity, y-direction
    vz = data[8]          # velocity, z-direction
    
    # Extract data for given star number
    list_range = len(time)
    x_star = np.zeros(list_range)
    y_star = np.zeros(list_range)
    z_star = np.zeros(list_range)
    vx_star = np.zeros(list_range)
    vy_star = np.zeros(list_range)
    vz_star = np.zeros(list_range)
    
    for i in range(len(time)):
            x_star[i] = x[i*total_stars + star_number]
            y_star[i] = y[i*total_stars + star_number]
            z_star[i] = z[i*total_stars + star_number]
            vx_star[i] = vx[i*total_stars + star_number]
            vy_star[i] = vy[i*total_stars + star_number]
            vz_star[i] = vz[i*total_stars + star_number]
                
    return time,m,x_star,y_star,z_star,vx_star,vy_star,vz_star

    

def plot_time(total_stars,star_number,time_step,subset):
    # Function that plots the results from VV and RK4 as a function of time, 
    # compared with the analytical solution.
    
    # Get data: Positions and velocities
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (total_stars,time_step)
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_%d_%.3f.txt' % (total_stars,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.3f.txt' % (total_stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.3f.txt' % (total_stars,time_step)
    
    tVV,mVV,xVV,yVV,zVV,vxVV,vyVV,vzVV = read_file(filename_verlet,star_number,total_stars)
    tRK4,mRK4,xRK4,yRK4,zRK4,vxRK4,vyRK4,vzRK4 = read_file(filename_RK4,star_number,total_stars)
    
    # Get data: Energies
    dataVV = np.loadtxt('../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_energy_%d_%.3f.txt' % (total_stars,time_step),unpack=True)
    dataRK4 = np.loadtxt('../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_energy_%d_%.3f.txt' % (total_stars,time_step),unpack=True)
    time = dataVV[0]
    EkVV = dataVV[1]
    EpVV = dataVV[2]
    EtotVV = dataVV[3]
    EkRK4 = dataRK4[1]
    EpRK4 = dataRK4[2]
    EtotRK4 = dataRK4[3]
    
    if subset == True:
        # Plot position: VV
        plt.figure()
        plt.title('Position, time step = %.3f' % time_step,size=12)
        plt.plot(tVV,xVV,label='VV')
        plt.plot(tVV[::5],xVV[::5],label='VV (subset 1:5)')
        plt.plot(tVV[::8],xVV[::8],label='VV (subset 1:10)')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$x$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
        
        # Plot position: RK4
        plt.figure()
        plt.title('Position, time step = %.3f' % time_step,size=12)
        plt.plot(tRK4,xRK4,label='RK4')
        plt.plot(tRK4[::5],xRK4[::5],label='RK4 (subset 1:5)')
        plt.plot(tRK4[::10],xRK4[::10],label='RK4 (subset 1:10)')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$x$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
        
        # Plot velocity: VV
        plt.figure()
        plt.title('Velocity, time step = %.3f' % time_step,size=12)
        plt.plot(tVV,vxVV,label='VV')
        plt.plot(tVV[::5],vxVV[::5],label='VV (subset 1:5)')
        plt.plot(tVV[::10],vxVV[::10],label='VV (subset 1:10)')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$v_x$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
        
        # Plot velocity: RK4
        plt.figure()
        plt.title('Velocity, time step = %.3f' % time_step,size=12)
        plt.plot(tRK4,vxRK4,label='RK4')
        plt.plot(tRK4[::5],vxRK4[::5],label='RK4 (subset 1:5)')
        plt.plot(tRK4[::10],vxRK4[::10],label='RK4 (subset 1:10)')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$v_x$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
        
        # Plot energy: VV
        plt.figure()
        plt.title('Total energy (Velocity-Verlet), time step = %.3f' % time_step,size=12)
        plt.plot(time,EkVV,'r',label='K')
        plt.plot(time[::5],EkVV[::5],'b',label='K (1:5)')
        plt.plot(time[::10],EkVV[::10],'g',label='K (1:10)')
        
        plt.plot(time,EpVV,'r',label='P')
        plt.plot(time[::5],EpVV[::5],'b',label='P (1:5)')
        plt.plot(time[::10],EpVV[::10],'g',label='P (1:10)')
        
        plt.plot(time,EtotVV,'r',label='Tot')        
        plt.plot(time[::5],EtotVV[::5],'b',label='Tot (1:5)')
        plt.plot(time[::10],EtotVV[::10],'g',label='Tot (1:10)')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$E$',size=14)
        plt.legend(loc=4,prop={'size':12})
        plt.show()
        
        # Plot energy: RK4
        plt.figure()
        plt.title('Total energy (RK4), time step = %.3f' % time_step,size=12)
        plt.plot(time,EkRK4,'r',label='K')
        plt.plot(time[::5],EkRK4[::5],'b',label='K (1:5)')
        plt.plot(time[::10],EkRK4[::10],'g',label='K (1:10)')
        
        plt.plot(time,EpRK4,'r',label='P')
        plt.plot(time[::5],EpRK4[::5],'b',label='P (1:5)')
        plt.plot(time[::10],EpRK4[::10],'g',label='P (1:10)')
        
        plt.plot(time,EtotRK4,'r',label='Tot')        
        plt.plot(time[::5],EtotRK4[::5],'b',label='Tot (1:5)')
        plt.plot(time[::10],EtotRK4[::10],'g',label='Tot (1:10)')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$E$',size=14)
        plt.legend(loc=4,prop={'size':12})
        plt.show()
        
    else:
        # Plot: Position
        plt.figure()
        plt.title('Position, time step = %.3f' % time_step,size=12)
        plt.plot(tVV,xVV,label='VV')
        plt.plot(tRK4,xRK4,label='RK4')
        plt.legend(loc=1,prop={'size':12})
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$x$',size=14)
        
        # Plot: Velocity
        plt.figure()
        plt.title('Velocity, time step = %.3f' % time_step,size=12)
        plt.plot(tVV,vxVV,label='Verlet')
        plt.plot(tRK4,vxRK4,label='RK4')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$v_x$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
        
        # Plot: Radial position
        rVV = np.sqrt(xVV**2 + yVV**2 + zVV**2)
        rRK4 = np.sqrt(xRK4**2 + yRK4**2 + zRK4**2)
    
        plt.figure()
        plt.title('Radial position, time step = %.2f' % time_step,size=12)
        plt.plot(tVV,rVV,label='Verlet')
        plt.plot(tRK4,rRK4,label='RK4')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$r$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
        
        # Plot: Energy
        plt.figure()
        plt.title('Total energy, time step = %.3f' % time_step,size=12)
        plt.plot(time,EkVV,'b--',label='VV kinetic')
        plt.plot(time,EpVV,'b.',label='VV potential')
        plt.plot(time,EtotVV,'b',label='VV total')
    
        plt.plot(time,EkRK4,'g--',label='RK4 kinetic')
        plt.plot(time,EpRK4,'g.',label='RK4 potential')
        plt.plot(time,EtotRK4,'g',label='RK4 total')
        plt.xlabel(r'$t$',size=14)
        plt.ylabel(r'$E$',size=14)
        plt.legend(loc=1,prop={'size':12})
        plt.show()
    
    return


def plot_orbits_VV(total_stars,time_step):
    
    # Get data
    filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_VV_%d_%.3f.txt' % (total_stars,time_step)
    #filename_verlet = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_VV_%d_%.3f.txt' % (total_stars,time_step)
    
    # Set up plot
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # Plot data
    for i in range(total_stars):
        t,m,x,y,z,vx,vy,vz = read_file(filename_verlet,i,total_stars)
        ax.plot(x,y,z, label='Star %d, VV' % i)
    
    radius = 20
    plt.xlim([-radius,radius])
    plt.ylim([-radius,radius])
    ax.set_zlim(-radius, radius)
    #ax.legend()
    plt.show()

    return
    

def plot_orbits_RK4(total_stars,time_step):
    
    # Get data
    filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_clang_64bit-Debug/cluster_RK4_%d_%.3f.txt' % (total_stars,time_step)
    #filename_RK4 = '../build-Project5-Desktop_Qt_5_5_0_MinGW_32bit-Debug/cluster_RK4_%d_%.3f.txt' % (total_stars,time_step)
    
    # Set up plot: RK4
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # Plot data
    for i in range(total_stars):
        t,m,x,y,z,vx,vy,vz = read_file(filename_RK4,i,total_stars)
        ax.plot(x,y,z, label='Star %d, RK4' % i)
    
    radius = 50
    #plt.xlim([-radius,radius])
    #plt.ylim([-radius,radius])
    #ax.legend()
    plt.show()

    return
    
    
def main(argv):
    total_stars = 100
    time_step = 0.001
    integration_points = 1000

    # Plot orbits
    #plot_orbits_VV(total_stars,time_step)
    
    # Plot results as a function of time
    #star_number = 0 # Plot for one particular star
    #subset = True # True = Print results as a function of time for different time steps, False = one time step
    #plot_time(total_stars,star_number,time_step,subset)
    
	
if __name__ == "__main__":
    main(sys.argv[1:]) 