#include "galaxy.h"
#include "star.h"
#include <iostream>
#include <cmath>
#include "time.h"

galaxy::galaxy()
{
    total_stars = 0;
    radius = 100;
    total_mass = 0;
}

galaxy::galaxy(double radi)
{
    total_stars = 0;
    radius = radi;
    total_mass = 0;
}

galaxy::galaxy(double radi,double mass)
{
    total_stars = 0;
    radius = radi;
    total_mass = mass;
}

void galaxy::add(star newstar)
{
    total_stars += 1;
    total_mass += newstar.mass;
    all_stars.push_back(newstar);
}

void galaxy::addM(star newstar)
{
    total_stars +=1;
    all_stars.push_back(newstar);
}

void galaxy::print_position(std::ofstream &output, int dimension, double time,int number)
{   // Writes mass, position and velocity to a file "output"
    if(dimension > 3 || dimension <= 0) dimension = 3;
    else{
        for(int i=0;i<number;i++){
            star &current = all_stars[i];
            output << time << "\t" << i+1 << "\t" << current.mass;
            for(int j=0;j<dimension;j++) output << "\t" << current.position[j];
            for(int j=0;j<dimension;j++) output << "\t" << current.velocity[j];
            output << std::endl;
        }
    }
}

void galaxy::RungeKutta4(int dimension, int integration_points, double final_time, bool stellar,int print_number)
{   // 4th order Runge-Kutta solver for a star cluster / spherical galaxy

    // Define time step
    double time_step = final_time/((double) integration_points);
    double time = 0.0;

    // Create file for data storage
    char *filename = new char[1000];
    if(stellar) sprintf(filename, "cluster_RK4_%d_%.2f.txt",total_stars,time_step); // If N-body cluster
    else sprintf(filename, "analytic_RK4_%d_%.2f.txt",integration_points,time_step); // If 1D 2-body analytic case
    std::ofstream output_file(filename);

    // Set up arrays
    double k1_x[total_stars][dimension],k2_x[total_stars][dimension],k3_x[total_stars][dimension],k4_x[total_stars][dimension];
    double k1_v[total_stars][dimension],k2_v[total_stars][dimension],k3_v[total_stars][dimension],k4_v[total_stars][dimension];

    double relative_position[3];
    double Fx,Fy,Fz;

    // Write initial values to file
    print_position(output_file,dimension,time,print_number);

    // Set up clock to measure the time usage
    clock_t start_RK4,finish_RK4;
    start_RK4 = clock();

    // START CALCULATIONS
    // Loop over time
    time += time_step;
    while(time<final_time){

        // k1
        for(int nr1=0;nr1<total_stars;nr1++){
            star &current = all_stars[nr1];
            Fx = Fy = Fz = 0.; // Reset forces for each run
            if(stellar){
                for(int nr2=nr1+1;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    for(int j=0;j<dimension;j++) relative_position[j] = -(other.position[j]-current.position[j]);
                    GravitationalForce_RK(relative_position[0],relative_position[1],relative_position[2],Fx,Fy,Fz,current.mass,other.mass);
                }
                for(int j=0;j<dimension;j++) k1_x[nr1][j] = time_step*current.velocity[j];
                k1_v[nr1][0] = time_step*Fx/current.mass;
                k1_v[nr1][1] = time_step*Fy/current.mass;
                k1_v[nr1][2] = time_step*Fz/current.mass;
            }
            else{
                Fx = -current.position[0];
                k1_x[nr1][0] = time_step*current.velocity[0];
                k1_v[nr1][0] = time_step*Fx/current.mass;
            }
        }

        // k2
        for(int nr1=0;nr1<total_stars;nr1++){
            star &current = all_stars[nr1];
            Fx = Fy = Fz = 0;
            if(stellar){
                for(int nr2=nr1+1;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    for(int j=0;j<dimension;j++) relative_position[j] = -((other.position[j]+k1_x[nr2][j]/2.)-(current.position[j]+k1_x[nr1][j]/2.));
                    GravitationalForce_RK(relative_position[0],relative_position[1],relative_position[2],Fx,Fy,Fz,current.mass,other.mass);
                }
                for(int j=0;j<dimension;j++) k2_x[nr1][j] = time_step*(current.velocity[j]+k1_v[nr1][j]/2.);
                k2_v[nr1][0] = time_step*Fx/current.mass;
                k2_v[nr1][1] = time_step*Fy/current.mass;
                k2_v[nr1][2] = time_step*Fz/current.mass;
            }
            else{
                Fx = -(current.position[0]+k1_x[nr1][0]/2.);
                k2_x[nr1][0] = time_step*(current.velocity[0]+k1_v[nr1][0]/2.);
                k2_v[nr1][0] = time_step*Fx/current.mass;
            }
        }

        // k3
        for(int nr1=0;nr1<total_stars;nr1++){
            star &current = all_stars[nr1];
            Fx = Fy = Fz = 0;
            if(stellar){
                for(int nr2=nr1+1;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    for(int j=0;j<dimension;j++) relative_position[j] = -((other.position[j]+k2_x[nr2][j]/2.)-(current.position[j]+k2_x[nr1][j]/2.));
                    GravitationalForce_RK(relative_position[0],relative_position[1],relative_position[2],Fx,Fy,Fz,current.mass,other.mass);
                }
                for(int j=0;j<dimension;j++) k3_x[nr1][j] = time_step*(current.velocity[j]+k2_v[nr1][j]/2.);
                k3_v[nr1][0] = time_step*Fx/current.mass;
                k3_v[nr1][1] = time_step*Fy/current.mass;
                k3_v[nr1][2] = time_step*Fz/current.mass;
            }
            else{
                Fx = -(current.position[0]+k2_x[nr1][0]/2.);
                k3_x[nr1][0] = time_step*(current.velocity[0]+k2_v[nr1][0]/2.);
                k3_v[nr1][0] = time_step*Fx/current.mass;
            }
        }

        // k4
        for(int nr1=0;nr1<total_stars;nr1++){
            star &current = all_stars[nr1];
            Fx = Fy = Fz = 0;
            if(stellar){
                for(int nr2=nr1+1;nr2<total_stars;nr2++){
                    star &other = all_stars[nr2];
                    for(int j=0;j<dimension;j++) relative_position[j] = -((other.position[j]+k3_x[nr2][j])-(current.position[j]+k3_x[nr1][j]));
                    GravitationalForce_RK(relative_position[0],relative_position[1],relative_position[2],Fx,Fy,Fz,current.mass,other.mass);
                }
                for(int j=0;j<dimension;j++) k4_x[nr1][j] = time_step*(current.velocity[j]+k3_v[nr1][j]);
                k4_v[nr1][0] = time_step*Fx/current.mass;
                k4_v[nr1][1] = time_step*Fy/current.mass;
                k4_v[nr1][2] = time_step*Fz/current.mass;
            }
            else{
                Fx = -(current.position[0]+k3_x[nr1][0]);
                k4_x[nr1][0] = time_step*(current.velocity[0]+k3_v[nr1][0]);
                k4_v[nr1][0] = time_step*Fx/current.mass;
            }
        }

        // Update position and velocity
        for(int nr1=0;nr1<total_stars;nr1++){
            star &current = all_stars[nr1];
            for(int j=0;j<dimension;j++){
                current.position[j] += (k1_x[nr1][j] + 2.*(k2_x[nr1][j] + k3_x[nr1][j]) + k4_x[nr1][j])/6.;
                current.velocity[j] += (k1_v[nr1][j] + 2.*(k2_v[nr1][j] + k3_v[nr1][j]) + k4_v[nr1][j])/6.;
            }
        }

        // Write current values to file and increase time
        print_position(output_file,dimension,time,print_number);
        time += time_step;
    }
    // Stop clock and print out time usage
    finish_RK4 = clock();
    std::cout << "Total time = " << "\t" << ((float)(finish_RK4 - start_RK4)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
    std::cout << "One time step = " << "\t" << ((float)(finish_RK4 - start_RK4)/CLOCKS_PER_SEC)/integration_points << " seconds" << std::endl; // print elapsed time

    // Close files
    output_file.close();
}

void galaxy::VelocityVerlet(int dimension, int integration_points, double final_time, bool stellar, int print_number)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    // Define time step
    double time_step = final_time/((double) integration_points);
    double time = 0.0;

    // Create files for data storage
    char *filename = new char[1000];
    if(stellar) sprintf(filename, "cluster_VV_%d_%.2f.txt",total_stars,time_step); // If N-body cluster
    else sprintf(filename, "analytic_VV_%d_%.2f.txt",integration_points,time_step); // If 1D 2-body analytic case
    std::ofstream output_file(filename);

    // Set up arrays
    double **acceleration = setup_matrix(total_stars,3);
    double **acceleration_new = setup_matrix(total_stars,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension

    // Write initial values to file
    print_position(output_file,dimension,time,print_number);

    // Set up clock to measure the time usage
    clock_t start_VV,finish_VV;
    start_VV = clock();

    // START CALCULATIONS
    // Loop over time
    time += time_step;
    while(time < final_time){

        // Loop over all stars
        for(int nr1=0; nr1<total_stars; nr1++){
            star &current = all_stars[nr1]; // Current star we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Calculate forces in each dimension
            if(stellar){
                for(int nr2=nr1+1; nr2<total_stars; nr2++){
                    star &other = all_stars[nr2];
                    GravitationalForce(current,other,Fx,Fy,Fz);
                }
            }
            else Fx = -current.position[0];

            // Acceleration in each dimension for current star
            acceleration[nr1][0] = Fx/current.mass;
            acceleration[nr1][1] = Fy/current.mass;
            acceleration[nr1][2] = Fz/current.mass;

            // Calculate new position for current star
            for(int j=0; j<dimension; j++) {
                current.position[j] += current.velocity[j]*time_step + 0.5*time_step*time_step*acceleration[nr1][j];
            }

            // Loop over all other stars
            if(stellar){
                for(int nr2=nr1+1; nr2<total_stars; nr2++){
                    star &other = all_stars[nr2];
                    GravitationalForce(current,other,Fxnew,Fynew,Fznew);
                }
            }
            else Fxnew = -current.position[0];

            // Acceleration each dimension exerted for current star
            acceleration_new[nr1][0] = Fxnew/current.mass;
            acceleration_new[nr1][1] = Fynew/current.mass;
            acceleration_new[nr1][2] = Fznew/current.mass;

            // Calculate new velocity for current star
             for(int j=0; j<dimension; j++) current.velocity[j] += 0.5*time_step*(acceleration[nr1][j] + acceleration_new[nr1][j]);
        }

        // Write current values to file and increase time
        print_position(output_file,dimension,time,print_number);
        time += time_step;
    }
    // Stop clock and print out time usage
    finish_VV = clock();
    std::cout << "Total time = " << "\t" << ((float)(finish_VV - start_VV)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
    std::cout << "One time step = " << "\t" << ((float)(finish_VV - start_VV)/CLOCKS_PER_SEC)/integration_points << " seconds" << std::endl; // print elapsed time

    // Close files
    output_file.close();

    // Clear memory
    delete_matrix(acceleration);
    delete_matrix(acceleration_new);
}

double ** galaxy::setup_matrix(int height,int width)
{   // Function to set up a 2D array

    // Set up matrix
    double **matrix;
    matrix = new double*[height];

    // Allocate memory
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];

    // Set values to zero
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

void galaxy::delete_matrix(double **matrix)
{   // Function to deallocate memory of a 2D array

    for (int i=0; i<total_stars; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void galaxy::GravitationalForce(star &current,star &other,double &Fx,double &Fy,double &Fz)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current star and all other stars
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);

    // Calculate the forces in each direction    
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/(r*r*r);
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/(r*r*r);
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/(r*r*r);
}

void galaxy::GravitationalForce_RK(double x_rel,double y_rel,double z_rel,double &Fx,double &Fy,double &Fz,double mass1,double mass2)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current star and all other stars
    double r = sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel);

    // Calculate the forces in each direction
    Fx -= this->G*mass1*mass2*x_rel/(r*r*r);
    Fy -= this->G*mass1*mass2*y_rel/(r*r*r);
    Fz -= this->G*mass1*mass2*z_rel/(r*r*r);
}

double galaxy::KineticEnergySystem()
{
    double KE = 0;
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr];
        KE += Current.KineticEnergy();
    }
    return KE;
}

double galaxy::PotentialEnergySystem()
{
    double PE = 0;
    for(int nr1=0;nr1<total_stars;nr1++){
        star &Current = all_stars[nr1];
        for(int nr2=nr1+1;nr2<total_stars;nr2++){
            star &Other = all_stars[nr2];
            PE += Current.PotentialEnergy(Other);
        }
    }
    return PE;
}

/*
bool galaxy::EnergyConservation(double previousEnergy)
{
    double TotalEnergy = this->KineticEnergySystem() + this->PotentialEnergySysyem();
    if((TotalEnergy - previousEnergy) < 0.05) return true;
    else return false;
}
*/

bool galaxy::Bound(star OneStar)
{
    //if(fabs(OneStar.position[0])>this->radius || fabs(OneStar.position[1])>this->radius || fabs(OneStar.position[2])>this->radius)
    //    return false;

    return ((OneStar.KineticEnergy() + this->PotentialEnergySystem()) < 0);
}

/*

void galaxy::Remove(int index)
{
    all_stars.erase(index);
    total_stars -= 1;
}

double galaxy::EnergyLoss()
{
    bool bound;
    double EnergyLoss = 0;
    int InitialStars = total_stars;
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr]
        bound = this->Bound(Current);
        if(!bound){
            this->Remove(nr);
            EnergyLoss -= Current.KineticEnergy();
        }
    }
    return EnergyLoss;
}

void galaxy::MergeStars()
{
    for(int nr1=0;nr1<total_stars;nr1++){
        star &Current = all_stars[nr1];
        for(int nr2=0;nr2<total_stars;nr2++){
            star &Other = all_stars[nr2];
            if(Current.distance(Other) < 0.001){
                Current.mass = Current.mass + Other.mass;
                this->Remove(nr2);
            }
        }
    }

}

*/
