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

