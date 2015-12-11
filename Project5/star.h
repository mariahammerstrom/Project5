#ifndef STAR_H
#define STAR_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class star
{
public:
    // Properties
    double mass;
    double position[3];
    double velocity[3];
    double potential;
    double kinetic;

    // Initializers
    star();
    star(double mas,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    double distance(star otherStar);
    double GravitationalForce(star otherStar, double Gconst);
    double Acceleration(star otherStar, double Gconst);
    //void merge(star star1,star star2);
    double KineticEnergy();
    double PotentialEnergy(star &otherStar, double Gconst, double epsilon);

};

#endif // STAR_H
