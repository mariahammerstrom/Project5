#ifndef GALAXY_H
#define GALAXY_H
#include "star.h"
#include <vector>
using std::vector;

class galaxy
{
public:
    int total_stars = 0; // initialize
    galaxy();
    vector<star> all_stars;
    void add(star newstar);
};

#endif // GALAXY_H
