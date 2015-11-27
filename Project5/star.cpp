#include "star.h"

star::star()
{
    mass = 1;
    position[0] = 1;
    position[1] = 0;
    position[2] = 0;
    velocity[0] = 0;
    velocity[1] = 0;
    velocity[2] = 0;
}

star::star(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
}


double star::distance(star otherStar)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = otherStar.position[0];
    y2 = otherStar.position[1];
    z2 = otherStar.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
}

double star::GravitationalForce(star otherStar)
{
    double r = this->distance(otherStar);
    if(r!=0) return G*this->mass*otherStar.mass/(r*r);
    else return 0;
}

double star::Acceleration(star otherStar)
{
    double r = this->distance(otherStar);
    if(r!=0) return this->GravitationalForce(otherStar)/this->mass/r;
    else return 0;
}

/*
void star::merge(star &star2)
{
    double m,vx,vy,vz;
    m = this->mass+star2.mass;
    vx = this->velocity[0]-star2.velocity[0];
    vy = this->velocity[1]-star2.velocity[1];
    vz = this->velocity[2]-star2.velocity[2];

    this->mass = m;
    this->velocity[0] = vx;
    this->velocity[1] = vy;
    this->velocity[2] = vz;

    star2.mass = 0;
    star2.position[0] = 0;
    star2.position[1] = 0;
    star2.position[2] = 0;
    star2.velocity[0] = 0;
    star2.velocity[1] = 0;
    star2.velocity[2] = 0;
}
*/
