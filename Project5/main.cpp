#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
//#include <mpi.h>
#include "star.h"
#include "galaxy.h"
#include "solver.h"

using namespace std;

void RK4(int dimension, int N, double final_time, double *x, double *v, double mass1, double mass2, ofstream &file);
void VV(int dimension, int N, double final_time, double *x_initial, double *v_initial, double mass1, double mass2, ofstream &file);
void print_initial(int dimension,double time_step, double final_time,double *x_initial,double *v_initial);
void print_final(int dimension, double *x_final, double *v_final);
void randomUniformSphere(double R0,double &x,double &y,double &z, default_random_engine *generator);
void Gaussian_distribution(double mean,double stddev,double &mass, default_random_engine *generator);

int main()
{
    double Msun,Mearth;
    Msun = 1.;
    Mearth = 1./332946;

    int N = 40;                     // No. of integration points
    double final_time = 10;         // End time of calculation
    double h = final_time/N;        // Time step

    // Set up initial conditions for each dimension
    int dimension = 3;
    double x_RK4[dimension],v_RK4[dimension],x_VV[dimension],v_VV[dimension];

    for(int i=0; i<dimension; i++){
        v_RK4[i] = v_VV[i] = 0.;
        x_RK4[i] = x_VV[i]= 0.3;
    }

    /*
    // Run algorithms
    if (dimension==1){
        cout << dimension << " DIMENSION:" << endl << endl;

        // Print out set up for the calculation
        print_initial(dimension,h,final_time,x_RK4,v_RK4);


        // 4TH-ORDER RUNGE-KUTTA
        cout << "RK4" << endl;

        // Set up file
        char filename_analytic_RK4[100];
        sprintf(filename_analytic_RK4, "RK4_analytic_%.3f.txt", h);
        ofstream file_analytic_RK4(filename_analytic_RK4);

        // Run calculation
        //for(double t=0; t<final_time; t+=h)
        RK4(dimension,N,final_time,x_RK4,v_RK4,Mearth,Msun,file_analytic_RK4);

        // Print out final result of the calculation
        print_final(dimension,x_RK4,v_RK4);

        // Close file
        file_analytic_RK4.close();


        // VELOCITY-VERLET
        cout << endl << "Velocity-Verlet" << endl;

        // Set up file
        char filename_analytic_VV[100];
        sprintf(filename_analytic_VV, "VV_analytic_%.3f.txt", h);
        ofstream file_analytic_VV(filename_analytic_VV);

        // Run calculation //for(double t=0; t<final_time; t+=h)
        VV(dimension,N,h,x_VV,v_VV,Mearth,Msun,file_analytic_VV);

        // Print out final result of the calculation
        print_final(dimension,x_VV,v_VV);

        // Close file
        file_analytic_VV.close();
    }

    else{
        cout << dimension << " DIMENSION:" << endl << endl;

        // Print out set up for the calculation
        print_initial(dimension,h,final_time,x_RK4,v_RK4);

        // 4TH-ORDER RUNGE-KUTTA
        cout << endl << "RK4" << endl;

        // Set up file
        char filename_RK4[100];
        sprintf(filename_RK4, "RK4_%.3f.txt", h);
        ofstream file_RK4(filename_RK4);

        // Run calculation //for(double t=0; t<final_time; t+=h)
        RK4(dimension,N,final_time,x_RK4,v_RK4,Mearth,Msun,file_RK4);

         // Print out final result of the calculation
        print_final(dimension,x_RK4,v_RK4);

        // Close file
        file_RK4.close();

        // VELOCITY-VERLET
        cout << endl << "Velocity-Verlet" << endl;

        // Set up file
        char filename_VV[100];
        sprintf(filename_VV, "VV_%.3f.txt", h);
        ofstream file_VV(filename_VV);

        // Run calculation
        //for(double t=0; t<final_time; t+=h)
        VV(dimension,N,h,x_VV,v_VV,Mearth,Msun,file_VV);

        // Print out final result of the calculation
        print_final(dimension,x_VV,v_VV);

        // Close file
        file_VV.close();
    }
    */

    // GALAXY MODEL
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    double R0 = 20; // Radius of galaxy
    int objects = 100; // Number of stars to be added in galaxy
    double m,x,y,z; // randomly distributed mass and position
    m = x = y = z = 0.0;
    double mean = 10.;
    double deviation = 1.;
    galaxy MM15(R0);
    for(int i=0;i<objects;i++){
        Gaussian_distribution(mean,deviation,m,&generator);
        cout << i << " " << m << endl;
        randomUniformSphere(R0,x,y,z,&generator);
        star stari(m,x,y,z,0,0,0);
        MM15.add(stari);
    }
    cout << "MM15 contains " << MM15.total_stars << " star(s)." << endl;

    char filename_galaxy[100];
    sprintf(filename_galaxy, "galaxy_%d_%.1f.txt",objects,R0);
    ofstream file_galaxy(filename_galaxy);
    MM15.print_position(file_galaxy,MM15.all_stars,dimension); // need vector<star> stars
    file_galaxy.close();

    return 0;
}

void RK4(int dimension,int N, double final_time, double *x, double *v, double mass1, double mass2,ofstream &file){
    // 4th order Runge-Kutta solver for two coupeled first ordered ODE

    double h = final_time/((double) N);

    double k1[dimension],k2[dimension],k3[dimension],k4[dimension]; // position
    double l1[dimension],l2[dimension],l3[dimension],l4[dimension]; // velocity

    double distance,r3,t;
    double xi = 0;
    double G = 4*M_PI*M_PI;


    for(int i=0;i<N;i++){
        t = i*h;
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

void VV(int dimension,int N,double final_time,double *x_initial,double *v_initial, double mass1, double mass2, ofstream &file){
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

    for(int i=1;i<N;i++){
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

void print_initial(int dimension,double time_step, double final_time,double *x_initial,double *v_initial){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << time_step << "; final time = " << final_time << endl;

    cout << "Initial position = ";
    for(int j=0;j<dimension;j++) cout << x_initial[j] << " ";
    cout << endl;

    cout << "Initial velocity = ";
    for(int j=0;j<dimension;j++) cout << v_initial[j] << " ";
    cout << endl;
}

void print_final(int dimension,double *x_final,double *v_final){
    // A function that prints out the final results of the calculation

    cout << "Final position = ";
    for(int j=0; j<dimension; j++) cout << x_final[j] << " ";
    cout << endl;

    cout << "Final velocity = ";
    for(int j=0; j<dimension; j++) cout << v_final[j] << " ";
    cout << endl;
}

void randomUniformSphere(double R0,double &x,double &y,double &z,default_random_engine *generator){
    // Random uniform number distribution that returns coordinates (x,y,z) in a sphere of radius R0.

    // Set up the uniform number generator
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //default_random_engine generator(seed);
    uniform_real_distribution<double> uniform_sphere(0.0,1.0);

    // Spherical coordinates
    double phi = 2*M_PI*uniform_sphere(*generator);
    double theta = acos(1 -2*uniform_sphere(*generator));
    double r = R0*pow(uniform_sphere(*generator),1./3);

    // Convert to cartesian coordinates
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);

}

void Gaussian_distribution(double mean,double stddev,double &mass, default_random_engine *generator){
    // Gaussian random number distribution that returns a mass around a mean with a given standard deviation.

    // Set up the uniform number generator
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //default_random_engine generator(seed);
    normal_distribution<double> normal_dist(mean,stddev);

    // Generate the mass
    mass = normal_dist(*generator);
    //cout << seed << endl;
}
