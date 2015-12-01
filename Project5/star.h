#ifndef STAR_H
#define STAR_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class star
{
public:
    // Astrophysical constants in units of solar masses, years, and lightyears
    double G = 4*M_PI*M_PI; // in AU
    //double G = 1.5608*1e-13; // in units of lightyears,years, and solar masses

    // Properties
    double mass;
    double position[3];
    double velocity[3];

    // Initializers
    star();
    star(double mas,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    double distance(star otherStar);
    double GravitationalForce(star otherStar);
    double Acceleration(star otherStar);
    //void merge(star star1,star star2);

};

#endif // STAR_H
