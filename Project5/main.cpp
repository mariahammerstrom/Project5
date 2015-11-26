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

void print_initial(int dimension,double time_step, double final_time,double *x_initial,double *v_initial);
void print_final(int dimension, double *x_final, double *v_final);
void randomUniformSphere(double R0,double &x,double &y,double &z, default_random_engine *generator);
void Gaussian_distribution(double mean,double stddev,double &mass, default_random_engine *generator);

int main()
{
    int N = 1000;               // No. of integration points
    double final_time = 10;     // End time of calculation
    int dimension;              // dimension of system
    bool stellar;               // true if star cluster


    // Testing code
    /*
    double Msun,Mearth;
    Msun = 1.;
    Mearth = 1./332946;

    double h = final_time/((double) N); // Time step

    stellar = false;
    dimension = 1;
    solver test;
    star star1(0.05,1,0,0,0,0,0);
    star star2(1,0,0,0,0,0,0);

    char filename_analytic_RK4[100];
    sprintf(filename_analytic_RK4, "RK4_analytic_%.3f.txt", h);
    ofstream file_analytic_RK4(filename_analytic_RK4);

    cout << "RK: " << endl;
    print_initial(dimension,h,final_time,star1.position,star1.velocity);
    test.RK4(dimension,N,final_time,star1,star2,file_analytic_RK4,stellar);
    print_final(dimension,star1.position,star1.velocity);
    cout << endl;
    file_analytic_RK4.close();

    star star3(0.05,1,0,0,0,0,0);
    star star4(1,0,0,0,0,0,0);

    char filename_analytic_VV[100];
    sprintf(filename_analytic_VV, "VV_analytic_%.3f.txt", h);
    ofstream file_analytic_VV(filename_analytic_VV);

    cout << "VV: " << endl;
    print_initial(dimension,h,final_time,star3.position,star3.velocity);
    test.VV(dimension,N,final_time,star3,star4,file_analytic_VV,stellar);
    print_final(dimension,star1.position,star3.velocity);
    file_analytic_VV.close();
    */


    // GALAXY (STAR CLUSTER) MODEL
    dimension = 3; stellar = true;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    double R0 = 20;         // Radius of galaxy
    int objects = 100;      // Number of stars to be added in galaxy
    double m,x,y,z;         // randomly distributed mass and position
    m = x = y = z = 0.0;

    // mean and standard deviation of stellar mass
    double mean = 10.;
    double deviation = 1.;

    galaxy MM15(R0);

    for(int i=0;i<objects;i++){
        Gaussian_distribution(mean,deviation,m,&generator);
        randomUniformSphere(R0,x,y,z,&generator);
        star stari(m,x,y,z,0,0,0);
        MM15.add(stari);
    }
    cout << "MM15 contains " << MM15.total_stars << " star(s)." << endl;

    // run system through RK4/VV, all data is written to file as we go
    MM15.RungeKutta4(dimension,N,final_time);
    // MM15.VelocityVerlet(dimension,N,final_time);

    return 0;
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
    normal_distribution<double> normal_dist(mean,stddev);

    // Generate the mass
    mass = normal_dist(*generator);
    //cout << seed << endl;
}
