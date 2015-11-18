#include "star.h"

star::star()
{

}

star::star(double mas, double radi, double x, double y, double z, double vx, double vy, double vz)
{
    mass = mas;
    radius = radi;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
}

double star::distance(star object1, star object2)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = object1.position[0];
    y1 = object1.position[1];
    z1 = object1.position[2];

    x2 = object2.position[0];
    y2 = object2.position[1];
    z2 = object2.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
}

double star::GravitationalForce(star star1, star star2)
{
    double force[3];
    double r = distance(star1,star2);
    for(i=0;i<3;i++)
        force[i] = -G*star1.mass*star2.mass/(r*r);
}

