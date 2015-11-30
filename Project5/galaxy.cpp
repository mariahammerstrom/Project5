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

void galaxy::RungeKutta4(int dimension,int N, double final_time, bool stellar)
{   // 4th order Runge-Kutta solver for a star cluster / spherical galaxy

    // Define time step
    double time_step = final_time/((double) N);
    double time = 0.;

    // Create file for data storage
    char *filename = new char[1000];
    if(stellar) sprintf(filename, "cluster_RK4_%d_%.2f.txt",total_stars,time_step); // If N-body cluster
    else sprintf(filename, "analytic_RK4_%d_%.2f.txt",N,time_step); // If 1D 2-body analytic case
    std::ofstream output_file(filename);

    // Set up arrays
    double **acceleration = setup_matrix(total_stars,dimension);
    double x1[total_stars][dimension],x2[total_stars][dimension],x3[total_stars][dimension],x4[total_stars][dimension];
    double v1[total_stars][dimension],v2[total_stars][dimension],v3[total_stars][dimension],v4[total_stars][dimension];
    double a1[total_stars][dimension],a2[total_stars][dimension],a3[total_stars][dimension],a4[total_stars][dimension];

    double dvdt,relative_position;
    double Fx,Fy,Fz;
    double position_temp_other[3];

    std::cout << "Time" << "\t" << "Star" << "\t" << "Mass"<< "\t" << "x" << "\t" << "y" << "\t" << "z" << "\t" << "vx" << "\t" << "vy" << "\t" << "vz" << "\t" << "Fx" << "\t" << "Fy" << "\t" << "Fz" << std::endl;
    std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;


    // Write initial values to file
    print_position(output_file,dimension,time);

    // Start evolving
    time = time_step;
    while(time<final_time){
        // k1
        for(int nr1=0;nr1<total_stars;nr1++){
            star &current = all_stars[nr1];
            dvdt = 0;
            Fx = Fy = Fz = 0.;
            for(int j=0;j<dimension;j++){
                if(stellar){
                    for(int nr2=0;nr2<total_stars;nr2++){
                        if(nr2!=nr1){
                            star &other = all_stars[nr2];
                            relative_position = (other.position[j]-current.position[j]);
                            dvdt += current.Acceleration(other)*relative_position;
                        }
                    }
                }else{
                    dvdt += current.position[j];
                }
                // x1[nr1][j] = current.position[j];
                v1[nr1][j] = time_step*current.velocity[j];
                a1[nr1][j] = -time_step*dvdt;
            }
        }
        // k2
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                if(stellar){
                    for(int nr2=0;nr2<total_stars;nr2++){
                        if(nr2!=nr){
                            star &other = all_stars[nr2];
                            relative_position = (other.position[j]+v1[nr2][j]/2.)-(current.position[j]+v1[nr][j]/2.);
                            dvdt += current.Acceleration(other)*relative_position;
                        }
                    }
                }else{
                    dvdt += current.position[j]+v1[nr][j]/2.;
                }
                v2[nr][j] = time_step*(current.velocity[j]+a1[nr][j]/2.);
                a2[nr][j] = -time_step*dvdt;
            }
        }
        // k3
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                if(stellar){
                    for(int nr2=0;nr2<total_stars;nr2++){
                        if(nr2!=nr){
                            star &other = all_stars[nr2];
                            relative_position = (other.position[j]+v2[nr2][j]/2.)-(current.position[j]+v2[nr][j]/2.);
                            dvdt += current.Acceleration(other)*relative_position;
                        }
                    }
                }else{
                    dvdt += current.position[j]+v2[nr][j]/2.;
                }
                v3[nr][j] = time_step*(current.velocity[j]+a2[nr][j]/2.);
                a3[nr][j] = -time_step*dvdt;
            }
        }
        // k4
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            dvdt = 0;
            for(int j=0;j<dimension;j++){
                if(stellar){
                    for(int nr2=0;nr2<total_stars;nr2++){
                        if(nr2!=nr){
                            star &other = all_stars[nr2];
                            relative_position = (other.position[j]+v3[nr2][j])-(current.position[j]+v3[nr][j]);
                            dvdt += current.Acceleration(other)*relative_position;
                        }
                    }
                }else{
                    dvdt += current.position[j]+v3[nr][j];
                }
                v4[nr][j] = time_step*(current.velocity[j]+a3[nr][j]);
                a4[nr][j] = -time_step*dvdt;
            }
        }

        // Update position and velocity
        for(int nr=0;nr<total_stars;nr++){
            star &current = all_stars[nr];
            for(int j=0;j<dimension;j++){
                current.position[j] += (v1[nr][j] + 2.*(v2[nr][j] + v3[nr][j]) + v4[j][nr])/6.;
                current.velocity[j] += (a1[nr][j] + 2.*(a2[nr][j] + a3[nr][j]) + a4[j][nr])/6.;
            }
        }

        // Write current values to file
        print_position(output_file,dimension,time);
        time += time_step;
    }
    output_file.close();
}

void galaxy::VelocityVerlet(int dimension,int N, double final_time, bool stellar)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    // Define time step
    double time_step = final_time/((double) N);
    double time = 0.0;

    // Create file for data storage
    char *filename = new char[1000];
    if(stellar) sprintf(filename, "cluster_VV_%d_%.2f.txt",total_stars,time_step); // If N-body cluster
    else sprintf(filename, "analytic_VV_%d_%.2f.txt",N,time_step); // If 1D 2-body analytic case
    std::ofstream output_file(filename);

    // Set up arrays
    double **acceleration = setup_matrix(total_stars,3);
    double **acceleration_new = setup_matrix(total_stars,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension

    std::cout << "Time" << "\t" << "Star" << "\t" << "Mass"<< "\t" << "x" << "\t" << "y" << "\t" << "z" << "\t" << "vx" << "\t" << "vy" << "\t" << "vz" << "\t" << "Fx" << "\t" << "Fy" << "\t" << "Fz" << std::endl;
    std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;


    // START CALCULATIONS
    // Loop over time
    while(time < final_time){

        // Loop over all stars
        for(int nr1=0; nr1<total_stars; nr1++){
            star &current = all_stars[nr1]; // Current star we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Loop over all other stars to calculate the force contributions from them
            for(int nr2=nr1+1; nr2<total_stars; nr2++){
                star &other = all_stars[nr2];

                // Calculate forces in each dimension
                if(stellar) GravitationalForce(current,other,Fx,Fy,Fz);
                else SpringForce(current,other,Fx);
            }
            // Acceleration in each dimension for current star
            acceleration[nr1][0] = Fx/current.mass;
            acceleration[nr1][1] = Fy/current.mass;
            acceleration[nr1][2] = Fz/current.mass;

            // Write to file
            if(dimension == 1 && nr1 == 0) output_file << time << "\t" << current.position[0] << "\t" << current.velocity[0] << "\t" << Fx << std::endl; // Analytical 1D 2-body case
            if(dimension != 1) print_position(output_file,dimension,time);

            // Calculate new position for current star
            for(int j=0; j<3; j++) current.position[j] = current.position[j] + current.velocity[j]*time_step + 0.5*time_step*time_step*acceleration[nr1][j];

            // Loop over all other stars
            for(int nr2=nr1+1; nr2<total_stars; nr2++){
                star &other = all_stars[nr2];

                // Calculate forces in each dimension using the updated position
                if(stellar) GravitationalForce(current,other,Fxnew,Fynew,Fznew);
                else SpringForce(current,other,Fxnew);
            }
            // Acceleration each dimension exerted for current star
            acceleration_new[nr1][0] = Fxnew/current.mass;
            acceleration_new[nr1][1] = Fynew/current.mass;
            acceleration_new[nr1][2] = Fznew/current.mass;

            // Calculate new velocity for current star
            for(int j=0; j<3; j++) current.velocity[j] = current.velocity[j] + 0.5*time_step*(acceleration[nr1][j] + acceleration_new[nr1][j]);

            // Print out values to inspect
            std::cout << time << "\t" << nr1+1 << "\t" << current.mass << "\t" << current.position[0] << "\t" << current.position[1] << "\t" << current.position[2] << "\t" << current.velocity[0] << "\t" << current.velocity[1] << "\t" << current.velocity[2] << "\t" << Fx << "\t" << Fy << "\t" << Fz << std::endl;
        }
        // Increase time step
        time += time_step;
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
    }
    // Close file
    output_file.close();

    // Clear memory
    delete_matrix(acceleration);
    delete_matrix(acceleration_new);
}

double ** galaxy::setup_matrix(int width,int height)
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

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j] - other.position[j];
    double r = sqrt(relative_distance[0]*relative_distance[0] + relative_distance[1]*relative_distance[1] + relative_distance[2]*relative_distance[2]);

    // Calculate the forces in each direction
    Fx += -G*current.mass*other.mass*relative_distance[0]/(r*r*r);
    Fy += -G*current.mass*other.mass*relative_distance[1]/(r*r*r);
    Fz += -G*current.mass*other.mass*relative_distance[2]/(r*r*r);
}

void galaxy::GravitationalForce_RK(double x1,double x2,double y1,double y2,double z1,double z2,double &Fx,double &Fy,double &Fz,double mass1,double mass2)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current star and all other stars
    double x_rel,y_rel,z_rel;

    for(int j = 0; j < 3; j++){
        x_rel = x1 - x2;
        y_rel = y1 - y2;
        z_rel = z1 - z2;
    }
    double r = sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel);

    // Calculate the forces in each direction
    Fx += -G*mass1*mass2*x_rel/(r*r*r);
    Fy += -G*mass1*mass2*y_rel/(r*r*r);
    Fz += -G*mass1*mass2*z_rel/(r*r*r);
}

void galaxy::SpringForce(star &current,star &other,double &Fx)
{   // Function that calculates the spring force between a light, moving object and a heavy, stationary object.

    // Calculate relative distance between current star and all other stars
    double relative_distance = current.position[0] - other.position[0];

    // Calculate the force
    Fx += -1.0*relative_distance;
}

void galaxy::SpringForce_RK(double x1,double x2,double &Fx)
{   // Function that calculates the spring force between a light, moving object and a heavy, stationary object.

    // Calculate relative distance between current star and all other stars
    double relative_distance = x1 - x2;

    // Calculate the force
    Fx += -1.0*relative_distance;
}
