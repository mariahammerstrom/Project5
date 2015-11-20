#include "galaxy.h"
#include "star.h"
#include <iostream>

galaxy::galaxy()
{
    total_stars = 0;
    radius = 100;
}

galaxy::galaxy(double radi)
{
    total_stars = 0;
    radius = radi;
}


void galaxy::add(star newstar)
{
    total_stars += 1;
    all_stars.push_back(newstar);
}

void galaxy::print_position(std::ofstream &output, vector<star> vec, int n, double time)
{   // writes mass, position and velocity to a file "output"
    if(n>3 || n<=0) n=3;
    for(int i=0;i<vec.size();i++){
        star &current = vec[i];
        output << time << "\t" << i+1 << "\t" << current.mass;
        for(int j=0;j<n;j++) output << "\t" << current.position[j];
        for(int j=0;j<n;j++) output << "\t" << current.velocity[j];
        output << std::endl;
    }
}

void galaxy::RungeKutta4(int dimension,int N, double t_crunch)
{
    char *filename = new char[1000];
    sprintf(filename, "cluster_RK4_%.2f.txt", t_crunch);
    std::ofstream output_file(filename);

    double Fg,relative_position;
    double k1[dimension],k2[dimension],k3[dimension],k4[dimension]; // position
    double l1[dimension],l2[dimension],l3[dimension],l4[dimension]; // velocity

    double time_step = t_crunch/((double) N);
    double time=0;
    print_position(output_file,this->all_stars,3,time);
    time = time_step;
    while(time<t_crunch){
        for(int nr=0;nr<total_stars;nr++){
            Fg = 0;
            star &current = all_stars[nr];
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    relative_position = (other.position[j]-current.position[j]);
                    Fg += current.GravitationalForce_r3(other)*relative_position;
                }
                k1[j] = time_step*current.velocity[j];
                l1[j] = -time_step*Fg;
            }
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    relative_position = ((other.position[j]+k1[j]/2.)-(current.position[j]+k1[j]/2.));
                    Fg += current.GravitationalForce_r3(other)*relative_position;
                }
                k2[j] = time_step*(current.velocity[j]+l1[j]/2.);
                l2[j] = -time_step*Fg;
            }
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    relative_position = ((other.position[j]+k2[j]/2.)-(current.position[j]+k2[j]/2.));
                    Fg += current.GravitationalForce_r3(other)*relative_position;
                }
                k3[j] = time_step*(current.velocity[j]+l2[j]/2.);
                l3[j] = -time_step*Fg;
            }
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    relative_position = ((other.position[j]+k3[j])-(current.position[j]+k3[j]));
                    Fg += current.GravitationalForce_r3(other)*relative_position;
                }
                k4[j] = time_step*(current.velocity[j]+l3[j]);
                l4[j] = -time_step*Fg;
            }
            for(int j=0;j<dimension;j++){
            current.position[j] += (k1[j] + 2.*(k2[j] + k3[j]) + k4[j])/6.;
            current.velocity[j] += (l1[j] + 2.*(l2[j] + l3[j]) + l4[j])/6.;
            }
        }
        print_position(output_file,this->all_stars,3,time);
        time += time_step;
    }
    output_file.close();
}

void galaxy::VelocityVerlet(int N, double t_crunch)
{
    for(int i=0;i<N;i++){
        this->all_stars;
    }
}

