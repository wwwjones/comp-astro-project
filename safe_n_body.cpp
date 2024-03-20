#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm> //for std::sort
#include <utility> //for std::pair
#include <chrono>
#include <Eigen/Core>
#include "particle.h"
#include "multipole.h"

using Matrix = Eigen::Matrix3d;
using Vector = Eigen::Vector3d;

unsigned int max_level;
unsigned int level;

void prep_data(std::ifstream &data, std::vector<double> &radii, std::vector<particle> &particles, double &M) {
    //Read in the data and create a vector of the radii
    double x, y, z, m, temp;
    while (!data.eof()){
    //for (int i = 0; i <= 10; ++i){
        particle* p = new particle;
        for (int j = 0; j < 10; ++j){
            if (j == 1) {
                data >> m;
                p->set_mass(m);
                M += m;
            }
            else if (j == 2){
                data >> x;
                p->set_x(x);
            }
            else if (j == 3){
                data >> y;
                p->set_y(y);
            }
            else if (j == 4){
                data >> z;
                p->set_z(z);
            }
            else data >> temp;
        }
        particles.push_back(*p);
        radii.push_back((p->position()).norm());
    }

    radii.pop_back();
    particles.pop_back();
    std::sort(radii.begin(), radii.end());
    std::cout << "radii size: " << radii.size() << "\nparticles size: " << particles.size() << std::endl;
}

void calc_distr(std::vector<double> radii) {
    double length = radii.size();

    //Could turn the counting into a function with a eq or exp bool
    //Count occurences in equidistant shells
    double r = 0.;
    unsigned int count = 0;
    double max = radii[length - 1];

    double bin_size = max / 10000.; //Number of bins
    std::vector<unsigned int> eq_sums;
    while (r < max) {
        r += bin_size;
        while (radii[count] < r && count < length) {
            ++count;
        }
        eq_sums.push_back(count);
    }

    //Count occurences in exponentially-spaced shells
    r = 0;
    count = 0;
    double x = -9.; //Experimentally determined start point (Could justify using differences between beginning data)
    std::vector<unsigned int> sums;
    std::vector<double> shells; //For some reason if I define this here I get a segmentation fault... writing over other stuff, but why?
    while (r < max) {
        r += std::exp(x);
        x += 0.01; //make much smaller for more measurements, 0.1?
        while (radii[count] < r && count < length) {
            ++count;
        }
        sums.push_back(count);
        shells.push_back(r);
        //compare with poisson distr. and calculate error
    }

    std::cout << "eq_sums total: " << eq_sums[eq_sums.size() - 1] << "\nsums total: " << sums[sums.size() - 1] << std::endl;

    //Print out a little graph (logarithmic)
    //for (unsigned int i = 5; i < sums.size(); i += 5){
    //    unsigned int stars = (sums[i] - sums[i - 5]) / 5;
    //    if (stars > 0) {
    //        for (unsigned int j = 0; j < stars; ++j)
    //            std::cout << '*';
    //        std::cout << std::endl;
    //    }
    //}

    //Export results for plotting
    std::ofstream results("results/results.txt");
    for (int i = 1; i < shells.size(); ++i){
        results << shells[i] << " ";
    }
    results << std::endl;
    for (int i = 1; i < sums.size(); ++i){
        results << sums[i] -  sums[i-1] << " ";
    }

    std::ofstream eq_results("results/eq_results.txt");
    for (int i = 1; i < eq_sums.size(); ++i){
        eq_results << eq_sums[i] -  eq_sums[i-1] << " ";
    }


    std::cout << "x: " << x << std::endl << "sums size: " << sums.size() << std::endl << "shells size: " << shells.size() << std::endl;
}


//Returns the index for the correct next cube pointer and updates the center
std::pair<int, Vector> pick_next_cube(Vector position, Vector center, double sidelength){
    //returns the id of the next cube
    //maybe a shorter implmentation? three checks to variables then a switch?
    std::pair<int, Vector> retvals;
    if(position.x() < center.x()){
        if(position.y() < center.y()){
            if(position.z() < center.z()){
                center.x() -= sidelength;
                center.y() -= sidelength;
                center.z() -= sidelength;
                retvals.first = 0;
                retvals.second = center;
                return retvals; //dlf
            }
            else {
                center.x() -= sidelength;
                center.y() -= sidelength;
                center.z() += sidelength;
                retvals.first = 1;
                retvals.second = center;
                return retvals; //ulf
            }
        }
        if(position.z() < center.z()){
            center.x() -= sidelength;
            center.y() += sidelength;
            center.z() -= sidelength;
            retvals.first = 2;
            retvals.second = center;
            return retvals; //dlb
        }
        else {
            center.x() -= sidelength;
            center.y() += sidelength;
            center.z() += sidelength;
            retvals.first = 3;
            retvals.second = center;
            return retvals; //uln
        }
    }
    if(position.y() < center.y()){
        if(position.z() < center.z()){
            center.x() += sidelength;
            center.y() -= sidelength;
            center.z() -= sidelength;
            retvals.first = 4;
            retvals.second = center;
            return retvals; //drf
        }
        else {
            center.x() += sidelength;
            center.y() -= sidelength;
            center.z() += sidelength;
            retvals.first = 5;
            retvals.second = center;
            return retvals; //urf
        }
    }
    if(position.z() < center.z()){
        center.x() += sidelength;
        center.y() += sidelength;
        center.z() -= sidelength;
        retvals.first = 6;
        retvals.second = center;
        return retvals; //drb
    }
    else {
        center.x() += sidelength;
        center.y() += sidelength;
        center.z() += sidelength;
        retvals.first = 7;
        retvals.second = center;
        return retvals; //urb
    }
}

void place(unsigned int i, std::vector<particle> &particles, Vector center, double sidelength, cube* current) {
    ++level;
    particle p = particles[i];
    int id = current->id();
    std::pair<int, Vector> idx_cntr; //for storing the index and the new center
    idx_cntr = pick_next_cube(p.position(), center, sidelength/4.);
    int index = idx_cntr.first;
    Vector next_center = idx_cntr.second;
    //Update this cubes info as we are passing through
    current->update_cube(p); //update the mass and CoM

    if (id >= 0){ //If the current cube has an occupant
        //Put our current particle in the right cube
        cube* next_cube = new cube(p, i, next_center, sidelength/2.);
        current->set_next(next_cube, index);
        //We have to re-place the inhabitant of this cube
        current->rmv_from_cube(p);//reverse the CoM and mass because it will be r>
        current->set_id(-1);
        --level;
        place(id, particles, center, sidelength, current);
    }
    else {
        cube* next = current->next(index);
        if (next){
            current = next;
            place(i, particles, next_center, sidelength/2., current);
        }
        else {
            cube* next_cube = new cube(p, i, next_center, sidelength/2.);
            current->set_next(next_cube, index);
            ++level;
        }
    }
}

//might actually be correct...
Matrix calculate_quadrupole(Vector com, cube* current){
    Matrix Q = Matrix::Zero();

    if (current->id() >= 0){
        for (int j = 0; j < 3; ++j){
            for (int i = 0; i < 3; ++i){
                //TODO: implement the right thing here
                //it's probably going to need a custom dirac-delta, easiest that way
                Q(i,j) = 1;
            }
        }
        return Q;
    }

    for (int i = 0; i < 8; ++i){
        cube* next = current->next(i);
        if (next) {
            Q += calculate_quadrupole(com, current->next(i));
        }
    }
    return Q;
}

//calculate it for current then
//there might be a very smart way to do this...
void calculate_properties(cube* current){
    Matrix Q = calculate_quadrupole(current->com(), current);
    current->set_quadrupole(Q);

    //radius calculation should be very similar, maybe even join into 1 function?
    double r = calculate_radius(current->com(), current);
    current->set_radius(r);

    for (int i = 0; i < 8; ++i){
        cube* next = current->next(i);
        if (next) calculate_properties(next);
    }
}

int main() {
    std::ifstream data ("data.txt");
    std::vector<double> radii;
    std::vector<particle> particles;
    double M = 0;

    //Read the data, create vectors of the radii and the masses and calculate total mass
    prep_data(data, radii, particles, M);

    //Calculate the distribution profile of the data and compare it to theory
    calc_distr(radii);

    //for (int i = 0; i < 10; ++i){
    //    std::cout << particles[i].x() << " " << particles[i].y() << " " << particles[i].z() << std::endl;
    //}

    //------------------------------------------------------------------------
    //--------------Calculate the N body forces of the particles--------------
    //------------------------------------------------------------------------

    std::vector<Vector> forces;
    //Make code after here into a function passing particles and forces
    Vector force, dif;
    double denom;
    double eps = 0.1; //softening - to be changed and better understood later
    double G = 1.; //gravitational constant, here == 1

    unsigned int its = particles.size();
    /*auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < its; ++i){
        if (i % 1000 == 0) std::cout << "Iteration " << i << " reached\n";
        force << 0, 0, 0;
        for (unsigned int j = 0; j < its; ++j){
            dif = particles[i].position() - particles[j].position();
            denom = std::pow((dif.squaredNorm() + std::pow(eps, 2)), 1.5);
            force += particles[i].mass() / denom * dif;
        }
        force *= -G;
        forces.push_back(force);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout << "Time for force calculation: " << time.count() << "s" << std::endl;

    std::ofstream results("results/forces.txt");

    for (int i = 0; i < its; ++i){
        results << particles[i].x() << " ";
    }
    results << std::endl;
    for (int i = 0; i < its; ++i){
        results << particles[i].y() << " ";
    }
    results << std::endl;
    for (int i = 0; i < its; ++i){
        results << particles[i].z() << " ";
    }
    results << std::endl;
    for (int i = 0; i < its; ++i){
        results << forces[i].norm() << " ";
    }

    */
    //for (int i = 0; i < 10; ++i){
    //    std::cout << forces[i].x() << " " << forces[i].y() << " " << forces[i].z() << " " << forces[i].norm() << std::endl;
    //}

    //------------------------------------------------------------------------
    //------------------------Build the multipole tree------------------------
    //------------------------------------------------------------------------

    //Place particles one at a time into the tree, all that should be needed is the position

    double sidelength = 1.;
    double max_radius = radii.back();
    while (sidelength < max_radius){
        sidelength *= 2.;
    }
    sidelength *= 2.;

    cube* head = new cube(particles[0], 0, Vector::Zero(), sidelength);

    //Build the tree
    for (unsigned int i = 1; i < its; ++i) {
        level = 0;
        place(i, particles, Vector::Zero(), sidelength, head);
        if (level > max_level) max_level = level;
    }

    //Calculate the quadrupoles
    calculate_properties(head);

    std::cout << "Particles placed: " << total_p << std::endl;
    std::cout << "Max sidelength: " << sidelength << std::endl;
    std::cout << "Number of levels: " << max_level << std::endl;
    std::cout << "Number of cubes: " << total_c << std::endl;
    std::cout << "Final Q: \n" << head->Q() << std::endl;

    return 0;
}
