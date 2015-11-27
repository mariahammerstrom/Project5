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

void galaxy::print_position(std::ofstream &output, int dimension, double time)
{   // writes mass, position and velocity to a file "output"
    if(dimension > 3 || dimension <= 0) dimension = 3;
    for(int i=0;i<total_stars;i++){
        star &current = all_stars[i];
        output << time << "\t" << i+1 << "\t" << current.mass;
        for(int j=0;j<dimension;j++) output << "\t" << current.position[j];
        for(int j=0;j<dimension;j++) output << "\t" << current.velocity[j];
        output << std::endl;
    }
}

void galaxy::RungeKutta4(int dimension,int N, double t_crunch)
{   // RK4 solver for a star cluster / spherical galaxy

    // if(stellar)
    char *filename = new char[1000];
    sprintf(filename, "cluster_RK4_%.2f.txt", t_crunch);
    std::ofstream output_file(filename);

    double dvdt,relative_position;
    double k1[total_stars][dimension],k2[total_stars][dimension],k3[total_stars][dimension],k4[total_stars][dimension];
    double l1[total_stars][dimension],l2[total_stars][dimension],l3[total_stars][dimension],l4[total_stars][dimension];

    double time_step = t_crunch/((double) N);
    double time = 0;

    // write initial values to file
    print_position(output_file,dimension,time);

    // start evolving
    time = time_step;
    while(time<t_crunch){
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    if(nr2!=nr){
                    star &other = all_stars[nr2];
                    relative_position = (other.position[j]-current.position[j]);
                    dvdt += current.Acceleration(other)*relative_position;}
                }
                k1[nr][j] = time_step*current.velocity[j];
                l1[nr][j] = -time_step*dvdt;
            }
        }
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    if(nr2!=nr){
                    star &other = all_stars[nr2];
                    relative_position = (other.position[j]+k1[nr2][j]/2.)-(current.position[j]+k1[nr][j]/2.);
                    dvdt += current.Acceleration(other)*relative_position;}
                }
                k2[nr][j] = time_step*(current.velocity[j]+l1[nr][j]/2.);
                l2[nr][j] = -time_step*dvdt;
            }
        }
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    if(nr2!=nr){
                    star &other = all_stars[nr2];
                    relative_position = ((other.position[j]+k2[nr][j]/2.)-(current.position[j]+k2[nr][j]/2.));
                    dvdt += current.Acceleration(other)*relative_position;}
                }
                k3[nr][j] = time_step*(current.velocity[j]+l2[nr][j]/2.);
                l3[nr][j] = -time_step*dvdt;
            }
        }
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                for(int nr2=0;nr2<total_stars;nr2++){
                    if(nr2!=nr){
                    star &other = all_stars[nr2];
                    relative_position = ((other.position[j]+k3[nr][j])-(current.position[j]+k3[nr][j]));
                    dvdt += current.Acceleration(other)*relative_position;}
                }
                k4[nr][j] = time_step*(current.velocity[j]+l3[nr][j]);
                l4[nr][j] = -time_step*dvdt;
            }
        }

            // update position and velocity
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            for(int j=0;j<dimension;j++){
            current.position[j] += (k1[nr][j] + 2.*(k2[nr][j] + k3[nr][j]) + k4[j][nr])/6.;
            current.velocity[j] += (l1[nr][j] + 2.*(l2[nr][j] + l3[nr][j]) + l4[j][nr])/6.;
            }
            //std::cout << nr << " " << current.position[0] << " " << current.position[1] << " " << current.position[2];
            //std::cout << " " << current.velocity[0]  << " " << current.velocity[1]  << " " << current.velocity[2];
            //std::cout << "k: " << k1[nr][0] << " " << k2[nr][0] << " " << k3[nr][0] << " " << k4[nr][0] << " ";
            //std::cout << "l: " << l1[nr][0] << " " << l2[nr][0] << " " << l3[nr][0] << " " << l4[nr][0] << " ";
            //std::cout << std::endl;
        }
        // std::cout << "-----------------" << std::endl;
        // write current values to file
        print_position(output_file,dimension,time);
        time += time_step;
    }
    output_file.close();
}

void galaxy::VelocityVerlet(int dimension,int N, double t_crunch)
{
    for(int i=0;i<N;i++){
        this->all_stars;
    }
}

