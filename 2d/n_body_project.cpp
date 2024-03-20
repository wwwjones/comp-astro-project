#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

// Read the position data and convert it into a vector of the radii
// Sort the list
// Find the max radius and map the interval onto [0, 1] so the center lambda is at 0.5
// Count the particles in bins of increasing radii (shells), determine the radii using a poisson distribution on the range [0, 1]
// Compare the results to the expected results according to the Hernquist paper

struct particle {
    double radius;
    double mass;
};

void prep_data(std::ifstream &data, std::vector<double> &radii, std::vector<double> &masses, double &M) {
    //Read in the data and create a vector of the radii
    double ri;
    double temp;
    while (!data.eof()){
        ri = 0.;
        for (int j = 0; j < 10; ++j){
            data >> temp;
            if (j == 1){
                masses.push_back(temp);
                M += temp;
            }
            if (j == 2 || j == 3 || j == 4){
                ri += temp * temp;
            }
        }
        ri = std::sqrt(ri);
        radii.push_back(ri);
    }

    radii.pop_back();
    masses.pop_back();
    //std::sort(radii.begin(), radii.end());
    std::cout << "radii size: " << radii.size() << "\nmasses size: " << masses.size() << std::endl;
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
    std::ifstream data ("../data.txt");
    std::vector<double> radii;
    std::vector<double> masses;
    double M = 0;

    //Read the data, create vectors of the radii and the masses and calculate total mass
    prep_data(data, radii, masses, M);

    std::vector<particle> particles;
    for (unsigned int i = 0; i < radii.size(); ++i){
        particle p = {radii[i], masses[i]};
        particles.push_back(p);
    }

    std::sort(particles.begin(), particles.end(), [](particle a, particle b) {return a.radius < b.radius;});

    //std::ofstream results("results/test.txt");
    //for (int i = 0; i < particles.size(); ++i){
    //    results << particles[i].radius << " " << particles[i].mass << std::endl;
    //}

    std::sort(radii.begin(), radii.end());

    //Calculate the distribution profile of the data and compare it to theory
    calc_distr(radii);

    //------------------------------------------------------------------------
    //--------------Calculate the N body forces of the particles--------------
    //------------------------------------------------------------------------

    std::vector<double> forces;
    double force, dif, denom;
    double eps = 0.1; //softening - to be changed and better understood later
    double G = 1.; //gravitational constant, here == 1
    //for convenience
    auto r = [&particles](int i){return particles[i].radius;};
    auto m = [&particles](int i){return particles[i].mass;};

    std::cout << "particles size: " << particles.size() << std::endl;

   auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < particles.size(); ++i){
        force = 0;
        for (unsigned int j = 0; j < particles.size(); ++j){
            dif = r(i) - r(j);
            denom = std::pow((std::pow(dif, 2) + std::pow(eps, 2)), 1.5);
            force += m(j) * dif / denom;
        }
        force *= -G;
        forces.push_back(force);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout << "Time for force calculation: " << time.count() << "s" << std::endl;

    std::ofstream results("results/forces.txt");
    for (int i = 0; i < forces.size(); ++i){
        results << r(i) << " ";
    }
    results << std::endl;
    for (int i = 0; i < forces.size(); ++i){
        results << forces[i] << " ";
    }

    return 0;
}
