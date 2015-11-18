#ifndef STAR_H
#define STAR_H


class star
{
public:
    // Astrophysical constants
    double G = 6.67259e-8; // cm3 g-1 s-2

    // Properties
    star();
    star(double mas,double radi,double x,double y,double z,double vx, double vy,double vz);
    double mass;
    double radius;
    double position[3];
    double velocity[3];
    double distance(star star1,star star2);
    double GravitationalForce(star star1, star star2);

};

#endif // STAR_H
