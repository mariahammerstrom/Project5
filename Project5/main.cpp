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

double dvdt(double x, double mass1, double mass2);
void RK4(int dimension, double dt, double final_time, double *x, double *v, double mass1, double mass2, ofstream &file);
void VV(int N, int dimension, double final_time, double *x_initial, double *v_initial, double mass1, double mass2, ofstream &file);

int main()
{
    double Msun,Mearth;
    Msun = 1.;
    Mearth = 1./332946;

    //double x = 1.; // initial position in AU
    //double v = 0.; // initial velocity in AU/yr
    int dimension = 3;
    double x[dimension],v[dimension];
    for(int i=0;i<dimension;i++){
        v[i] = 0.;
        x[i] = 0.;
    }
    x[0]=1.;

    double time_step,final_time;
    time_step = 0.01;
    final_time = 10;

    // Runge-Kutta 4th order
    ofstream file_RK("RK.txt");
    cout << "RK4" << endl;
    cout << "Time step = " << time_step << "; final time = " << final_time << endl;
    cout << "Initial position = ";
    for(int j=0;j<dimension;j++) cout << x[j] << " ";
    cout << endl;
    cout << "Initial velocity = ";
    for(int j=0;j<dimension;j++) cout << v[j] << " ";
    cout << endl;
    RK4(dimension,time_step,final_time,x,v,Mearth,Msun,file_RK);
    cout << "Final position = ";
    for(int j=0;j<dimension;j++) cout << x[j] << " ";
    cout << endl;
    cout << "Final velocity = ";
    for(int j=0;j<dimension;j++) cout << v[j] << " ";
    cout << endl;
    file_RK.close();


    // Velocity-Verlet
    for(int i=0;i<dimension;i++){
        v[i] = 0.;
        x[i] = 0.;
    }
    x[0]=1.;
    int N = final_time/time_step;
    double h = final_time/((double) N);
    ofstream file_VV("VV.txt");
    cout << endl << "Velocity-Verlet" << endl;
    cout << "Time step = " << h << "; final time = " << final_time << endl;
    cout << "Initial position = ";
    for(int j=0;j<dimension;j++) cout << x[j] << " ";
    cout << endl;
    cout << "Initial velocity = ";
    for(int j=0;j<dimension;j++) cout << v[j] << " ";
    cout << endl;
    VV(N,dimension,final_time,x,v,Mearth,Msun,file_VV);
    cout << "Final position = ";
    for(int j=0;j<dimension;j++) cout << x[j] << " ";
    cout << endl;
    cout << "Final velocity = ";
    for(int j=0;j<dimension;j++) cout << v[j] << " ";
    cout << endl;
    file_VV.close();

    return 0;
}

double dvdt(double x, double mass1, double mass2){
    double G = 4*M_PI*M_PI; // Gravitational constant
    double y = 0;
    double r = sqrt(x*x+y*y); // Distance between the two objects
    return -G*mass1*mass2*x/(r*r*r);
}

void RK4(int dimension,double h, double final_time, double *x, double *v, double mass1, double mass2,ofstream &file){
    // 4th order Runge-Kutta solver for two coupeled first ordered ODE

    double k1[dimension],k2[dimension],k3[dimension],k4[dimension]; // position
    double l1[dimension],l2[dimension],l3[dimension],l4[dimension]; // velocity

    double distance,r3;
    double xi = 0;
    double G = 4*M_PI*M_PI;


    for(double t=0;t<final_time;t+=h){
        xi = 0;
        for(int i=0;i<dimension;i++){
            xi += x[i]*x[i];
        }
        distance = sqrt(xi);
        r3 = distance*distance*distance;


        for(int j=0;j<dimension;j++){
            k1[j] = h*v[j];
            l1[j] = -h*G*mass1*mass2*x[j]/r3; // h*dvdt(x[j],mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k2[j] = h*(v[j]+l1[j]/2.);
            l2[j] = -h*G*mass1*mass2*(x[j]+k1[j]/2.)/r3; // h*dvdt(x[j]+k1[j]/2.,mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k3[j] = h*(v[j]+l2[j]/2.);
            l3[j] = -h*G*mass1*mass2*(x[j]+k2[j]/2.)/r3; // h*dvdt(x[j]+k2[j]/2.,mass1,mass2);
        }
        for(int j=0;j<dimension;j++){
            k4[j] = h*(v[j]+l3[j]);
            l4[j] = -h*G*mass1*mass2*(x[j]+k3[j])/r3; // h*dvdt(x[j]+k3[j],mass1,mass2);
        }

        // Update new position and velocity, write to file
        file << t << "\t";
        for(int j=0;j<dimension;j++){
            x[j] += (k1[j] + 2.*(k2[j] + k3[j]) + k4[j])/6.;
            v[j] += (l1[j] + 2.*(l2[j] + l3[j]) + l4[j])/6.;
            file << x[j] << "\t" << v[j] << "\t";
        }
        file << endl;
    }
}

void VV(int N,int dimension,double final_time,double *x_initial,double *v_initial, double mass1, double mass2, ofstream &file){
    // Velocity-Verlet

    double h = final_time/((double) N);
    double x[dimension*N];
    double v[dimension*N];
    double t;

    for(int j=0;j<dimension;j++){
        x[j] = x_initial[j];       // initial position
        v[j] = v_initial[j];       // initial velocity
        x[j+dimension] = x[j] + v[j]*h*h; // Euler
    }

    double distance,r3;
    double xi = 0;
    double G = 4*M_PI*M_PI;

    for(int i=1; i<N; i+=1){
        xi = 0;
        for(int i=0;i<dimension;i++){
            xi += x[i]*x[i];
        }
        distance = sqrt(xi);
        r3 = distance*distance*distance;

        t = i*h;
        for(int j=0;j<dimension;j++){
            x[j+(i+1)*dimension] = 2*x[j+i*dimension] - x[j+(i-1)*dimension] - h*h*G*mass1*mass2*x[j+i*dimension]/r3;
            v[j+i*dimension] = (x[j+(i+1)*dimension] - x[j+(i-1)*dimension])/(2*h);
        }

        // Write to file
        file << t << "\t";
        for(int j=0;j<dimension;j++)
            file << x[j+i*dimension] << "\t" << v[j+i*dimension] << "\t";
        file << endl;
    }
    for(int j=0;j<dimension;j++){
        x_initial[j] = x[dimension*N-(dimension-j)];       // final position
        v_initial[j] = v[dimension*N-(dimension-j)];       // final velocity
    }
}


