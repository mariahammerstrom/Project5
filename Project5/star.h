#ifndef STAR_H
#define STAR_H


class star
{
public:
    double mass;
    double radius;
    double position[3];
    double velocity[3];
    star();
    star(double mas,double radi,double x,double y,double z,double vx, double vy,double vz);

};

#endif // STAR_H
