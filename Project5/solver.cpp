#include "solver.h"
#include "star.h"
#include <cmath>

solver::solver()
{

}

void solver::RK4(int dimension,double h, double final_time,star &object1,star &object2){
    // 4th order Runge-Kutta solver for two coupeled first ordered ODE

    double k1[dimension],k2[dimension],k3[dimension],k4[dimension]; // position
    double l1[dimension],l2[dimension],l3[dimension],l4[dimension]; // velocity

    double r,Fg; // distance and gravitational force between the two objects

    for(double t=0;t<final_time;t+=h){
        r = distance(object1,object2);
        Fg = GravitationalForce(object1,object2);

        for(int j=0;j<dimension;j++){
            k1[j] = h*object1.velocity[j];
            l1[j] = -h*Fg*object1.position[j]/r; // h*dvdt(x[j],mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k2[j] = h*(object1.velocity[j]+l1[j]/2.);
            l2[j] = -h*Fg*(object1.position[j]+k1[j]/2.)/r; // h*dvdt(x[j]+k1[j]/2.,mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k3[j] = h*(object1.velocity[j]+l2[j]/2.);
            l3[j] = -h*Fg*(object1.position[j]+k2[j]/2.)/r; // h*dvdt(x[j]+k2[j]/2.,mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k4[j] = h*(object1.velocity[j]+l3[j]);
            l4[j] = -h*Fg*(object1.position[j]+k3[j])/r; // h*dvdt(x[j]+k3[j],mass1,mass2);
        }

        // Update new position and velocity
        for(int j=0;j<dimension;j++){
            object1.position[j] += (k1[j] + 2.*(k2[j] + k3[j]) + k4[j])/6.;
            object1.velocity[j] += (l1[j] + 2.*(l2[j] + l3[j]) + l4[j])/6.;
        }
    }
}


void solver::VV(int dimension,int N,double final_time,star &object1,star &object2){
    // Velocity-Verlet

    double h = final_time/((double) N);
    double x[dimension*N]; // position
    double v[dimension*N]; // velocity

    for(int j=0;j<dimension;j++){
        x[j] = object1.position[j];       // initial position
        v[j] = object1.velocity[j];       // initial velocity
        x[j+dimension] = x[j] + v[j]*h*h; // Euler
    }

    double r,Fg;

    for(int i=1; i<N; i+=1){
        r = distance(object1,object2);
        Fg = GravitationalForce(object1,object2);

        for(int j=0;j<dimension;j++){
            x[j+(i+1)*dimension] = 2*x[j+i*dimension] - x[j+(i-1)*dimension] - h*h*Fg*x[j+i*dimension]/r;
            v[j+i*dimension] = (x[j+(i+1)*dimension] - x[j+(i-1)*dimension])/(2*h);
        }
    }

    // Update position and velocity
    for(int j=0;j<dimension;j++){
        object1.position[j] = x[dimension*N-(dimension-j)];
        object1.velocity[j] = v[dimension*N-(dimension-j)];
    }
}
