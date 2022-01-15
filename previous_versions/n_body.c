#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <Eigen/Core>

// Read the position data and convert it into a vector of the radii
// Sort the list
// Find the max radius and map the interval onto [0, 1] so the center lambda is at 0.5
// Count the particles in bins of increasing radii (shells), determine the radii using a poisson distribution on the range [0, 1]
// Compare the results to the expected results according to the Hernquist paper

struct particle {
    Eigen::Vector3d position;
    double mass;
};

void prep_data(std::ifstream &data, std::vector<double> &radii, std::vector<particle> &particles, double &M) {
    //Read in the data and create a vector of the radii
    double ri;
    double temp;
    while (!data.eof()){
        particle* p = new particle;
        for (int j = 0; j < 10; ++j){
            if (j == 1) {
                data >> p->mass;
                M += p->mass;
            }
            else if (j == 2) data >> p->position[0];
            else if (j == 3 )data >> p->position[1];
            else if (j == 4 )data >> p->position[2];
            else data >> temp;
        }
        particles.push_back(*p);
        radii.push_back((p->position).norm());
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
    //    std::cout << particles[i].position.x() << " " << particles[i].position.y() << " " << particles[i].position.z() << std::endl;
    //}
    //------------------------------------------------------------------------
    //--------------Calculate the N body forces of the particles--------------
    //------------------------------------------------------------------------

    std::vector<Eigen::Vector3d> forces;
    Eigen::Vector3d force, dif;
    double denom;
    double eps = 0.1; //softening - to be changed and better understood later
    double G = 1.; //gravitational constant, here == 1
    //for convenience
    //auto r = [&particles](int i){return particles[i].radius;};
    //auto m = [&particles](int i){return particles[i].mass;};

    //std::cout << "particles size: " << particles.size() << std::endl;

    unsigned int its = particles.size();
    auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < its; ++i){
        force << 0, 0, 0;
        for (unsigned int j = 0; j < its; ++j){
            dif = particles[i].position - particles[j].position;
            denom = std::pow((dif.squaredNorm() + std::pow(eps, 2)), 1.5);
            force += particles[i].mass / denom * dif;
        }
        force *= -G;
        forces.push_back(force);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout << "Time for force calculation: " << time.count() << "s" << std::endl;

    std::ofstream results("results/forces.txt");

    for (int i = 0; i < its; ++i){
        results << particles[i].position.x() << " ";
    }
    results << std::endl;
    for (int i = 0; i < its; ++i){
        results << particles[i].position.y() << " ";
    }
    results << std::endl;
    for (int i = 0; i < its; ++i){
        results << particles[i].position.z() << " ";
    }
    results << std::endl;
    for (int i = 0; i < its; ++i){
        results << forces[i].norm() << " ";
    }


    //for (int i = 0; i < 10; ++i){
    //    std::cout << forces[i].x() << " " << forces[i].y() << " " << forces[i].z() << " " << forces[i].norm() << std::endl;
    //}

    return 0;
}
