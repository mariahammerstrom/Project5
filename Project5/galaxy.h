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
    double radius,total_mass,G;
    int total_stars;
    vector<star> all_stars;
    double totalKinetic;
    double totalPotential;

    // constants

    // initializers
    galaxy();
    galaxy(double radi);

    // functions
    void add(star newstar);
    void addM(star newstar);
    void GravitationalConstant();
    void print_position(std::ofstream &output, int dimension, double time, int number);
    void print_energy(std::ofstream &output, double time, double epsilon);
    void RungeKutta4(int dimension, int integration_points, double final_time, bool stellar, bool simple, int print_number, double epsilon);
    void VelocityVerlet(int dimension, int integration_points, double final_time, bool stellar, bool simple, int print_number, double epsilon);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
    void GravitationalForce(star &current, star &other, double &Fx, double &Fy, double &Fz, double epsilon);
    void GravitationalForce_RK(double x_rel, double y_rel, double z_rel, double &Fx, double &Fy, double &Fz, double mass1, double mass2);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    double EnergyLoss();
    bool Bound(star OneStar);

};

#endif // GALAXY_H
