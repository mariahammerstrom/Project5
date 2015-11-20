#ifndef STAR_H
#define STAR_H
#define _USE_MATH_DEFINES
#include <cmath>


class star
{
public:
    // Astrophysical constants in units of Msun, years, and lightyears
    double G = 4*M_PI*M_PI; // in AU

    // Initializers
    star();
    star(double mas,double x,double y,double z,double vx, double vy,double vz);

    // Properties
    double mass;
    double position[3];
    double velocity[3];

    // Functions
    double distance(star otherStar);
    double GravitationalForce(star otherStar);
    double GravitationalForce_r3(star otherStar);
    //void merge(star star1,star star2);

};

#endif // STAR_H
