#include "galaxy.h"
#include "star.h"

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

void galaxy::print_position(std::ofstream &output, vector<star> vec, int n)
{   // writes mass, position and velocity to a file "output"
    if(n>3 || n<=0) n=3;
    for(int i=0;i<vec.size();i++){
        star &current = vec[i];
        output << i+1 << "\t" << current.mass;
        for(int j=0;j<n;j++) output << "\t" << current.position[j];
        for(int j=0;j<n;j++) output << "\t" << current.velocity[j];
        output << std::endl;
    }
}

