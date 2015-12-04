#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "star.h"
#include "galaxy.h"

using namespace std;

void print_initial(int dimension, double time_step, double final_time, double *x_initial, double *v_initial, int N);
void print_final(int dimension, double *x_final, double *v_final);
void randomUniformSphere(double R0,double &x,double &y,double &z, default_random_engine *generator);
void Gaussian_distribution(double mean,double stddev,double &mass, default_random_engine *generator);

int main()
{

    // Which system would you like to run?
    // 1) Analytical test with box on a spring; set spring_test = true
    // 2) Two planets in a gravitational field; set binary = true
    // 3) Full star cluster in a gravitational field; set cluster = true

    bool spring_test,binary,cluster;

    spring_test = false;
    binary = false;
    cluster = true;

    int integration_points;  // No. of integration points
    double final_time;       // End time of calculation
    int dimension;           // No. of spatial dimensions

    bool force; // false = run program with analytical spring force, true = run program with gravitational potential


    // ANALYTIC: Spring force
    if(spring_test){
        dimension = 1;
        integration_points = 100;
        final_time = 100;
        force = false;

        cout << "ANALYTICAL" << endl;
        cout << "Time step: " << final_time/((double) integration_points) << ", integration points: " << integration_points << endl;

        // RK4 test
        star star1(1.,1.,0,0,0,0,0);
        galaxy testRK;
        testRK.add(star1);
        testRK.RungeKutta4(dimension,integration_points,final_time,force,1);

        // VV test
        star star2(1.,1.,0,0,0,0,0);
        galaxy testVV;
        testVV.add(star2);
        testVV.VelocityVerlet(dimension,integration_points,final_time,force,1);
    }

    // Binary stars
    if(binary){
        integration_points = 10000;
        final_time = 10000.;
        dimension = 3;
        double time_step = final_time/((double) integration_points);
        double x[3],v[3];
        force = true;

        cout << "BINARY SYSTEM" << endl;

        star star1(1.,60.,60.,60.,1.0,0.,0.);
        star star2(10.,0.,0.,0.,0.,0.,0.);

        galaxy binary_rk(100.0);
        binary_rk.add(star1);
        binary_rk.add(star2);

        for(int j=0;j<dimension;j++){
            x[j] = star1.position[j];
            v[j] = star1.velocity[j];
        }

        galaxy binary_vv(5.0);
        binary_vv.add(star1);
        binary_vv.add(star2);

        print_initial(dimension,time_step,final_time,x,v,integration_points);

        // Evolution of binary system
        cout << endl << "RK4: " << endl;
        binary_rk.RungeKutta4(dimension,integration_points,final_time,force,1);

        for(int j=0;j<dimension;j++){
            x[j] = binary_rk.all_stars[0].position[j];
            v[j] = binary_rk.all_stars[0].velocity[j];
        }
        print_final(dimension,x,v);

        cout << endl << "VV:" << endl;
        binary_vv.VelocityVerlet(dimension,integration_points,final_time,force,1);

        for(int j=0;j<dimension;j++){
            x[j] = binary_vv.all_stars[0].position[j];
            v[j] = binary_vv.all_stars[0].velocity[j];
        }
        print_final(dimension,x,v);

    }



    // GALAXY (STAR CLUSTER) MODEL
    if(cluster){
        dimension = 3;
        integration_points = 1000;
        final_time = 10;
        force = true;

        cout << "CLUSTER" << endl;
        cout << "Time step: " << final_time/((double) integration_points) << endl;
        cout << "Integration points: " << integration_points << endl;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);

        double R0 = 20.;         // Radius of galaxy
        int objects = 100;      // Number of stars to be added in galaxy
        double m,x,y,z;         // randomly distributed mass and position
        m = x = y = z = 0.0;

        // mean and standard deviation of stellar mass
        double mean = 10.;
        double deviation = 1.;

        galaxy MM15_rk(R0);
        galaxy MM15_vv(R0);

        for(int i=0;i<objects;i++){
            Gaussian_distribution(mean,deviation,m,&generator);
            randomUniformSphere(R0,x,y,z,&generator);
            star stari(m,x,y,z,0,0,0);
            MM15_rk.add(stari);
            MM15_vv.add(stari);
        }
        cout << "The star cluster MM15 contains " << MM15_rk.total_stars << " star(s)." << endl;

        // run system through RK4/VV, all data is written to file as we go
        MM15_rk.RungeKutta4(dimension,integration_points,final_time,force,5);
        MM15_vv.VelocityVerlet(dimension,integration_points,final_time,force,5);
    }

    return 0;
}



void print_initial(int dimension,double time_step, double final_time,double *x_initial,double *v_initial, int N){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << time_step << "; final time = " << final_time << "; integration points = " << N << endl;

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
}
