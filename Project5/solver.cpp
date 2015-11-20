#include "solver.h"
#include <cmath>
#include "star.h"

solver::solver()
{

}

void solver::RK4(int dimension,int N, double final_time,star &object1,star &object2){
    // 4th order Runge-Kutta solver for two coupeled ODE in "dimension" dimensions

    double k1[dimension],k2[dimension],k3[dimension],k4[dimension]; // position
    double l1[dimension],l2[dimension],l3[dimension],l4[dimension]; // velocity

    double Fg; // gravitational force between the two objects
    double h = final_time/((double) N); // time step

    for(int i=0;i<N;i++){
        Fg = object1.GravitationalForce_r3(object2);

        for(int j=0;j<dimension;j++){
            k1[j] = h*object1.velocity[j];
            l1[j] = -h*Fg*object1.position[j]; // h*dvdt(x[j],mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k2[j] = h*(object1.velocity[j]+l1[j]/2.);
            l2[j] = -h*Fg*(object1.position[j]+k1[j]/2.); // h*dvdt(x[j]+k1[j]/2.,mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k3[j] = h*(object1.velocity[j]+l2[j]/2.);
            l3[j] = -h*Fg*(object1.position[j]+k2[j]/2.); // h*dvdt(x[j]+k2[j]/2.,mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k4[j] = h*(object1.velocity[j]+l3[j]);
            l4[j] = -h*Fg*(object1.position[j]+k3[j]); // h*dvdt(x[j]+k3[j],mass1,mass2);
        }

        // Update new position and velocity
        for(int j=0;j<dimension;j++){
            object1.position[j] += (k1[j] + 2.*(k2[j] + k3[j]) + k4[j])/6.;
            object1.velocity[j] += (l1[j] + 2.*(l2[j] + l3[j]) + l4[j])/6.;
        }
    }
}


void solver::VV(int dimension,int N,double final_time,star &object1,star &object2){
    // Velocity-Verlet solver for two coupeled ODEs in "dimension" dimensions

    double h = final_time/((double) N); // time step
    double x[dimension*N];              // position
    double v[dimension*N];              // velocity

    for(int j=0;j<dimension;j++){
        x[j] = object1.position[j];       // initial position
        v[j] = object1.velocity[j];       // initial velocity
        x[j+dimension] = x[j] + v[j]*h*h; // Euler guess
    }

    double Fg; // gravitational force between the two objects

    for(int i=1; i<N; i+=1){
        Fg = object1.GravitationalForce_r3(object2);

        for(int j=0;j<dimension;j++){
            x[j+(i+1)*dimension] = 2*x[j+i*dimension] - x[j+(i-1)*dimension] - h*h*Fg*x[j+i*dimension];
            v[j+i*dimension] = (x[j+(i+1)*dimension] - x[j+(i-1)*dimension])/(2*h);
        }
    }

    // Update new position and velocity
    for(int j=0;j<dimension;j++){
        object1.position[j] = x[dimension*N-(dimension-j)];
        object1.velocity[j] = v[dimension*N-(dimension-j)];
    }
}
