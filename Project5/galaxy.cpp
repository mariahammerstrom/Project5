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
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
}

galaxy::galaxy(double radi)
{
    total_stars = 0;
    radius = radi;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
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

void galaxy::GravitationalConstant()
{
    G = (4*M_PI*M_PI/32)*radius*radius*radius/total_mass;
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

void galaxy::print_energy(std::ofstream &output, double time,double epsilon)
{   // Writes energies to a file "output"

    this->KineticEnergySystem();
    this->PotentialEnergySystem(epsilon);
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr];
        output << time << "\t" << nr << "\t";
        output << Current.kinetic << "\t" << Current.potential << std::endl;
    }
}

void galaxy::RungeKutta4(int dimension, int integration_points, double final_time, bool stellar, bool simple, int print_number,double epsilon)
{   // 4th order Runge-Kutta solver for a star cluster / spherical galaxy

    if(!simple) GravitationalConstant();

    // Define time step
    double time_step = final_time/((double) integration_points);
    double time = 0.0;
    double loss = 0.; // Possible energy loss

    // Create file for data storage
    char *filename = new char[1000];
    char *filenameE = new char[1000];
    if(stellar){
        sprintf(filename, "cluster_RK4_%d_%.3f.txt",total_stars,time_step); // If N-body cluster
        sprintf(filenameE, "cluster_RK4_energy_%d_%.3f.txt",total_stars,time_step);
    }
    else sprintf(filename, "analytic_RK4_%d_%.3f.txt",integration_points,time_step); // If 1D 2-body analytic case
    std::ofstream output_file(filename);
    std::ofstream output_energy(filenameE);

    // Set up arrays
    double k1_x[total_stars][dimension],k2_x[total_stars][dimension],k3_x[total_stars][dimension],k4_x[total_stars][dimension];
    double k1_v[total_stars][dimension],k2_v[total_stars][dimension],k3_v[total_stars][dimension],k4_v[total_stars][dimension];

    double relative_position[3];
    double Fx,Fy,Fz;

    // Write initial values to file
    print_position(output_file,dimension,time,print_number);
    print_energy(output_energy,time,epsilon);

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

        print_energy(output_energy,time,epsilon);
        loss += EnergyLoss();
    }
    // Stop clock and print out time usage
    finish_RK4 = clock();
    std::cout << "Total time = " << "\t" << ((float)(finish_RK4 - start_RK4)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
    std::cout << "One time step = " << "\t" << ((float)(finish_RK4 - start_RK4)/CLOCKS_PER_SEC)/integration_points << " seconds" << std::endl; // print elapsed time

    // Close files
    output_file.close();
    output_energy.close();
}

void galaxy::VelocityVerlet(int dimension, int integration_points, double final_time, bool stellar, bool simple, int print_number, double epsilon)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    if(!simple) GravitationalConstant();
    //std::cout << "G = " << G << std::endl;

    // Define time step
    double time_step = final_time/((double) integration_points);
    double time = 0.0;
    double loss = 0.; // Possible energy loss
    int lostStars[integration_points];

    // Create files for data storage
    char *filename = new char[1000];
    char *filenameE = new char[1000];
    char *filenameB = new char[1000];
    char *filenameLost = new char[1000];
    if(stellar){
        sprintf(filename, "cluster_VV_%d_%.3f.txt",total_stars,time_step); // If N-body cluster
        sprintf(filenameE, "cluster_VV_energy_%d_%.3f.txt",total_stars,time_step);
        sprintf(filenameB,"cluster_bound_%d_%.3f.txt",total_stars,time_step);
        sprintf(filenameLost,"cluster_lost_%d_%.3f.txt",total_stars,time_step);
    }
    else sprintf(filename, "analytic_VV_%d_%.3f.txt",integration_points,time_step); // If 1D 2-body analytic case
    std::ofstream output_file(filename);
    std::ofstream output_energy(filenameE);
    std::ofstream output_bound(filenameB);
    std::ofstream output_lost(filenameLost);

    // Set up arrays
    double **acceleration = setup_matrix(total_stars,3);
    double **acceleration_new = setup_matrix(total_stars,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension

    // Write initial values to file
    print_position(output_file,dimension,time,print_number);
    print_energy(output_energy,time,epsilon);

    int n = 0;
    lostStars[n] = 0;
    output_lost << time << "\t" << lostStars[n] << std::endl;
    n+=1;

    // Set up clock to measure the time usage
    clock_t start_VV,finish_VV;
    start_VV = clock();

    // START CALCULATIONS
    // Loop over time
    time += time_step;
    while(time < final_time){
        lostStars[n] = 0;

        // Loop over all stars
        for(int nr1=0; nr1<total_stars; nr1++){
            star &current = all_stars[nr1]; // Current star we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Calculate forces in each dimension
            if(stellar){
                for(int nr2=nr1+1; nr2<total_stars; nr2++){
                    star &other = all_stars[nr2];
                    GravitationalForce(current,other,Fx,Fy,Fz,epsilon);
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
                    GravitationalForce(current,other,Fxnew,Fynew,Fznew,epsilon);
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

        // Energy conservation

        // Write current values to file and increase time
        print_position(output_file,dimension,time,print_number);
        print_energy(output_energy,time,epsilon);

        loss += EnergyLoss();

        for(int nr=0;nr<total_stars;nr++){
            star &Current = all_stars[nr];
            if(!(this->Bound(Current))){
                lostStars[n] += 1;
            }
        }
        output_lost << time << "\t" << lostStars[n] << std::endl;
        n += 1;
        time += time_step;
    }
    // Stop clock and print out time usage
    finish_VV = clock();
    std::cout << "Total time = " << "\t" << ((float)(finish_VV - start_VV)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
    std::cout << "One time step = " << "\t" << ((float)(finish_VV - start_VV)/CLOCKS_PER_SEC)/integration_points << " seconds" << std::endl; // print elapsed time

    //loss = EnergyLoss();
    std::cout << "Total energyloss due to unbound stars: " << loss << std::endl;

    double boundStars = 0;
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr];
        if(this->Bound(Current)){
            output_bound << nr << std::endl;
            boundStars += 1;
        }
    }
    std::cout << "There are " << boundStars << " bound stars at the end of the run" << std::endl;

    // Close files
    output_file.close();
    output_energy.close();
    output_bound.close();
    output_lost.close();

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

void galaxy::GravitationalForce(star &current,star &other,double &Fx,double &Fy,double &Fz,double epsilon)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current star and all other stars
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    double smoothing = epsilon*epsilon*epsilon;

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r) + smoothing);
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r) + smoothing);
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r) + smoothing);
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

void galaxy::KineticEnergySystem()
{
    totalKinetic = 0;
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr];
        Current.kinetic = Current.KineticEnergy();
    }
}

void galaxy::PotentialEnergySystem(double epsilon)
{
    totalPotential = 0;
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr];
        Current.potential = 0;
    }
    for(int nr1=0;nr1<total_stars;nr1++){
        star &Current = all_stars[nr1];
        for(int nr2=nr1+1;nr2<total_stars;nr2++){
            star &Other = all_stars[nr2];
            Current.potential += Current.PotentialEnergy(Other,G,epsilon);
            Other.potential += Other.PotentialEnergy(Current,G,epsilon);
        }
    }
}

bool galaxy::Bound(star OneStar)
{
    return ((OneStar.kinetic + OneStar.potential) < 0.0);
}


double galaxy::EnergyLoss()
{
    bool bound;
    vector<int> indices;
    double EnergyLoss = 0;
    for(int nr=0;nr<total_stars;nr++){
        star &Current = all_stars[nr];
        bound = this->Bound(Current);
        if(!bound){
            indices.push_back(nr);
        }
    }
    for(int i=0;i<indices.size();i++) EnergyLoss += all_stars[indices[i]].KineticEnergy();
    return EnergyLoss;
}

/*
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
