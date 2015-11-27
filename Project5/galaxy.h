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
    void RungeKutta4(int dimension, int N, double t_crunch, bool stellar);
    void VelocityVerlet(int dimension, int N, double t_crunch);
};

#endif // GALAXY_H
