#ifndef GALAXY_H
#define GALAXY_H
#include "star.h"
#include <vector>
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
};

#endif // GALAXY_H
