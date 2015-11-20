#ifndef SOLVER_H
#define SOLVER_H
#include "star.h"


class solver
{
public:
    solver();
    friend class star;
    void RK4(int dimension,int N, double final_time,star &object1,star &object2); // Runge-Kutta 4th order method
    void VV(int dimension,int N,double final_time,star &object1,star &object2); // Velocity-Verlet method
};

#endif // SOLVER_H
