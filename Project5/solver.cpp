#include "solver.h"
#include "star.h"
#include <cmath>

solver::solver()
{

}

/*
double solver::derivatives(star object1, star object2)
{
    double r,dv1dt[3],dv2dt[3],dx1dt[3],dx2dt[3];
    r = distance(object1,object2);
    for(int i=0;i<3;i++){
        dx1dt[i] = object1.velocity[i];
        dx2dt[i] = object2.velocity[i];
        dv1dt[i] = -G*object1.mass*object2.mass*object1.position[i]/(r*r*r);
        dv2dt[i] = -G*object1.mass*object2.mass*object2.position[i]/(r*r*r);
    }
    return dx1dt,dx2dt,dv1dt,dv2dt;
}

solver::dxdt(star object){
    double dxdt[3];
    for(int i=0;i<3;i++)
        dxdt[i] = object.velocity[i];
    return dxdt;
}

solver::dvdt(double mass,double x,double y,double z,double MassOther){
    double dvdt[3];
    double r = sqrt(x*x+y*y+z*z);
    dvdt[0] = -G*mass*MassOther*x/(r*r*r);
    dvdt[1] = -G*mass*MassOther*y/(r*r*r);
    dvdt[2] = -G*mass*MassOther*z/(r*r*r);
    return dvdt;
}


solver::RungeKutta4(double h,double final_time,double time_step,star &object1)
{
    double k1,k2,k3,k4;
    for(int time=0;time<final_time;time+=time_step){
        k1 = dvdt(mass,x,y,z,MassOther);
        k2 = dvdt(mass,x+k1/2.,y+k1/2.,z+k1/2.,MassOther);
        k3 = dvdt(mass,x+k2/2.,y+k2/2.,z+k2/2.,MassOther);
        k4 = dvdt(mass,x+k3.,y+k3.,z+k3.,MassOther);
        x += h*(k1 + 2.*(k2 + k3) + k4)/6.;
        y += h*(k1 + 2.*(k2 + k3) + k4)/6.;
        z += h*(k1 + 2.*(k2 + k3) + k4)/6.;
    }
}


solver::VelocityVerlet(double h,int N,double x, double vx) // not finished
{
    for(int i=1;i<N;i++){
        vx =
                x += h + vx;
    }
}
*/
