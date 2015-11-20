#include "star.h"
//#include <cmath>
//#include <vector>
//using std::vector;

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

/*
void star::manyStars(int number)
{
    vector<star> many;
    for(int i=0;i<number;i++){
        star newStar;
        many.push_back(newStar);
    }
}*/

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
    return -G*this->mass*otherStar.mass/(r*r);
}

double star::GravitationalForce_r3(star otherStar)
{
    double r = this->distance(otherStar);
    return -G*this->mass*otherStar.mass/(r*r*r);
}

/*
void star::merge(star star1, star star2)
{
    mass = star1.mass+star2.mass;
    star new_star(mass,x,y,z,vx,vy,vz);

}
*/

