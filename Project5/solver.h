#ifndef SOLVER_H
#define SOLVER_H
#include "star.h"
#include <fstream>


class solver
{
public:
    solver();

    friend class star;

    void RK4(int dimension, int N, double final_time, star &star1, star star2, std::ofstream &file, bool &stellar); // Runge-Kutta 4th order method
    void VV(int dimension, int N, double final_time, star &star1, star star2, std::ofstream &file, bool &stellar); // Velocity-Verlet method
};

#endif // SOLVER_H
