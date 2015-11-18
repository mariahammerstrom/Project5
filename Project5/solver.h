#ifndef SOLVER_H
#define SOLVER_H
#include "star.h"


class solver
{
public:
    solver();
    double distance(star object1,star object2);
    double TwoBody(double t, double x);
    double RungeKutta4(double h,int N,double time,double x); // Runge-Kutta 4th order method
    double VelocityVerlet(double t,int N); // Velocity-Verlet method
};

#endif // SOLVER_H
