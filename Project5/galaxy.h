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

    // initializers
    galaxy();
    galaxy(double radi);

    // properties
    double radius;
    int total_stars;
    vector<star> all_stars;

    // functions
    void add(star newstar);
    void print_position(std::ofstream &output, vector<star> vec, int n);
};

#endif // GALAXY_H
