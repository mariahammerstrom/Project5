#ifndef GALAXY_H
#define GALAXY_H
#include "star.h"
#include <vector>
#include <fstream>
using std::vector;

class galaxy
{
public:
    friend class star;

    // properties
    double radius,total_mass;
    int total_stars;
    vector<star> all_stars;

    // constants
    double G = (4*M_PI*M_PI/32)*radius*radius*radius/(total_mass);

    // initializers
    galaxy();
    galaxy(double radi);
    galaxy(double radi,double mass);

    // functions
    void add(star newstar);
    void addM(star newstar);
    void print_position(std::ofstream &output, int dimension, double time, int number);
    void RungeKutta4(int dimension, int integration_points, double final_time, bool stellar, int print_number);
    void VelocityVerlet(int dimension, int integration_points, double final_time, bool stellar,int print_number);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
    void GravitationalForce(star &current,star &other,double &Fx,double &Fy,double &Fz);
    void GravitationalForce_RK(double x_rel, double y_rel, double z_rel, double &Fx, double &Fy, double &Fz, double mass1, double mass2);
    double KineticEnergySystem();
    double PotentialEnergySystem();
    // bool EnergyConservation();
    bool Bound(star OneStar);
    // void RemoveStar();

};

#endif // GALAXY_H
