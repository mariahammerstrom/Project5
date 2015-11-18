#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
//#include <mpi.h>
//#include "star.h"
//#include "galaxy.h"

using namespace std;

double dvdt(double x, double v, double mass1, double mass2);
void RK4(double dt, double final_time, double &x, double &v, double mass1, double mass2);
void VV(double dt, double final_time, double &x, double &v, double mass1, double mass2);

int main()
{

    double Msun,Mearth;
    Msun = 1.;
    Mearth = 1./332946;

    double x = 1; // initial position in AU
    double v = 0; // initial velocity in AU/yr

    double time_step,final_time;

    time_step = 1;
    final_time = 10;

    // Runge-Kutta 4th order
    cout << "RK4" << endl;
    cout << "Time step = " << time_step << "; final time = " << final_time << endl;
    cout << "Initial position = " << x << "; initial velocity = " << v << endl;
    RK4(time_step,final_time,x,v,Mearth,Msun);
    cout << "Final position = " << x << "; final velocity = " << v << endl;

    // Velocity-Verlet
    x=1,v=0;
    cout << endl << "Velocity-Verlet" << endl;
    cout << "Time step = " << time_step << "; final time = " << final_time << endl;
    cout << "Initial position = " << x << "; initial velocity = " << v << endl;
    VV(time_step,final_time,x,v,Mearth,Msun);
    cout << "Final position = " << x << "; final velocity = " << v << endl;



    return 0;
}

double dvdt(double x, double v, double mass1, double mass2){
    double G = 4*M_PI*M_PI; // Gravitational constant
    double y = 0;
    double r = sqrt(x*x+y*y); // Distance between the two objects
    return -G*mass1*mass2*x/(r*r*r);
}

void RK4(double dt, double final_time, double &x, double &v, double mass1, double mass2){
    // 4th order Runge-Kutta solver for two coupeled first ordered ODE
    double k1,k2,k3,k4,l1,l2,l3,l4;
    for(int t=0;t<final_time;t+=dt){
        k1 = dt*v;
        l1 = dt*dvdt(x,v,mass1,mass2);
        k2 = dt*(v+l1/2.);
        l2 = dt*dvdt(x+k1/2.,v+l1/2.,mass1,mass2);
        k3 = dt*(v+l2/2.);
        l3 = dt*dvdt(x+k2/2.,v+l2/2.,mass1,mass2);
        k4 = dt*(v+l3);
        l4 = dt*dvdt(x+k3,v+l3,mass1,mass2);

        // Update new position and velocity
        x += (k1 + 2.*(k2 + k3) + k4)/6.;
        v += (l1 + 2.*(l2 + l3) + l4)/6.;
    }
}


void VV(double dt, double final_time, double &x,double &v, double mass1, double mass2){
    // Velocity-Verlet
    double xmin = 1;
    double temp;
    for(int t=0;t<final_time;t+=dt){
        temp = x;
        x = 2*x - xmin + dt*dt*dvdt(x,v,mass1,mass2);
        v = (x-xmin)/(2*dt);
        xmin = temp;
    }
}


