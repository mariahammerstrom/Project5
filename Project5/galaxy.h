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
    double radius;
    int total_stars;
    vector<star> all_stars;

    // initializers
    galaxy();
    galaxy(double radi);

    // functions
    void add(star newstar);
    void print_position(std::ofstream &output, int dimension, double time);
    void RungeKutta4(int dimension, int N, double final_time, bool stellar);
    void VelocityVerlet(int dimension, int N, double final_time, bool stellar);
    double **setup_matrix(int width,int height);
    void delete_matrix(double **matrix);
    void GravitationalForce(star &current,star &other,double &Fx,double &Fy,double &Fz);
    void GravitationalForce_RK(double x1,double x2,double y1,double y2,double z1,double z2,double &Fx,double &Fy,double &Fz,double mass1,double mass2);
    void SpringForce(star &current,star &other,double &Fx);
    void SpringForce_RK(double x1,double x2,double &Fx);

};

#endif // GALAXY_H
