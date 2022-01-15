#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

// Read the position data and convert it into a vector of the radii
// Sort the list
// Find the max radius and map the interval onto [0, 1] so the center lambda is at 0.5
// Count the particles in bins of increasing radii (shells), determine the radii using a poisson distribution on the range [0, 1]
// Compare the results to the expected results according to the Hernquist paper

int main() {
    std::ifstream data ("data.txt");
    std::vector<double> radii;
    double temp;
    double M = 0;

    //Read in the data and create a vector of the radii
    while (!data.eof()){
    //for (int i = 0; i < 10; ++i){
        double ri = 0.;
        for (int j = 0; j < 10; ++j){
            data >> temp;
            if (j == 1){
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
    unsigned int length = radii.size();
    std::sort(radii.begin(), radii.end());

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
    //std::vector<double> shelllllls; For some reason if I define this here I get a segmentation fault... writing over other stuff, but why?
    while (r < max) {
        r += std::exp(x);
        x += 0.01; //make much smaller for more measurements, 0.1?
        while (radii[count] < r && count < length) {
            ++count;
        }
        sums.push_back(count);
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

    //Save the shell radii
    std::vector<double> shells;
    r = 0;
    //x = -9.;
    for (double i = -9; i < x; i += 0.01){
        r += std::exp(i);
        shells.push_back(r);
    }
    //while (r < max) {
    //    r += std::exp(x);
    //    x += 0.1;
    //    shells.push_back(r);
    //}

    //Export results for plotting
    std::ofstream max_r("results/max_radius.txt");
    max_r << max;

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

    return 0;
}
