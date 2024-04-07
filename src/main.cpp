#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm> //for std::sort
#include <utility>   //for std::pair
#include <chrono>
#include <string>
#include <iomanip>
#include <omp.h>
#include <Eigen/Core>
#include "particle.hpp"
#include "cell.hpp"
#include "oct_tree.hpp"

// [L] = 4.90445e+17 m = 15.8941 parsecs
// [M] = 2e+28 kg = 0.01 solar masses
// [T] = 2.97286e+17 s = 9.42043e+09 years

using Matrix = Eigen::Matrix3d;
using Vector = Eigen::Vector3d;

const int NUM_PARTICLES = 100000;
int max_level;
int level;
double G = 1.;
std::vector<long double> units(4, 0.);

// determine the units used in this simulation, where we want the half-mass radius to be equal to 3 parsecs
void det_units()
{
    const long double G = 6.67408e-11;
    const long double solar_mass = 2e30;
    const long double parsec = 3.0857e16;
    const long double year = 60. * 60. * 24. * 365.25;
    const long double Rhm = 0.188749;
    const long double mass = 92.4259;
    long double L, M, T, F;
    long double time_u, mass_u, length_u;
    long double length;

    length = Rhm; // Want Rhm to be equal to 3 parsecs

    // new unit = real value * unit/ unitless value
    // mass * mass_u == (mass/100) * solar_mass
    mass_u = (mass / 100) * solar_mass / mass;
    // length * length_u == 3 * parsec
    length_u = 3 * parsec / length;
    // we now use G to determine the new time unit
    time_u = std::sqrt((1. / G) * std::pow(length_u, 3) / (mass_u));

    L = length_u / parsec;
    M = mass_u / solar_mass;
    T = time_u / year;
    F = length_u * mass_u / (time_u * time_u) / 1e9;

    std::cout << "[L] = " << length_u << " m = " << L << " parsecs" << std::endl;
    std::cout << "[M] = " << mass_u << " kg = " << M << " solar masses" << std::endl;
    std::cout << "[T] = " << time_u << " s = " << T << " years" << std::endl;
    std::cout << "[F] = " << F << " GN" << std::endl;

    units = {L, M, T, F};
}

// read in the data from an external file and create a vector of the particles as well as determine some other values of interest
void prep_data(std::ifstream &data, std::vector<Particle> &particles, double &M, double &quarter_mass_r, double &tcross)
{
    // Read in the data and create a vector of the radii
    double m, x, y, z, vx, vy, vz, temp;
    while (!data.eof())
    {
        Particle p;
        data >> temp >> m >> x >> y >> z >> vx >> vy >> vz >> temp >> temp;
        M += m;
        p.set_mass(m);
        p.set_x(x);
        p.set_y(y);
        p.set_z(z);
        p.set_vx(vx);
        p.set_vy(vy);
        p.set_vz(vz);
        p.set_radius(p.position().norm());
        particles.push_back(p);
    }
    // Always reads in an extra line of random data
    M -= m;
    particles.pop_back();

    double N = particles.size();
    /*
        Vector CoM = Vector::Zero();
        for (int i = 0; i < N; ++i){
            CoM += particles[i].mass() * particles[i].position();
        }
        CoM /= M;
        std::cout << CoM << std::endl;
        for (int i = 0; i < N; ++i){
            particles[i].set_radius((particles[i].position() - CoM).norm());
        }
    */
    std::sort(particles.begin(), particles.end(), [](Particle &a, Particle &b)
              { return (a.radius() < b.radius()); });

    double half_mass = 0;
    int idx = 0;
    while (half_mass < M / 4.)
    {
        half_mass += particles[idx++].mass();
    }
    quarter_mass_r = particles[idx - 1].radius();

    while (half_mass < M / 2.)
    {
        half_mass += particles[idx++].mass();
    }
    double Rhm = particles[idx - 1].radius();

    double vc = std::sqrt(G * half_mass / Rhm);

    tcross = Rhm / vc;

    double trelax = N * tcross / (8. * std::log(N));

    double mips = std::cbrt(4. / 3. * M_PI * Rhm * Rhm * Rhm / (N / 2));

    std::cout << "Total mass: " << M << "\nQuarter-mass radius: " << quarter_mass_r << "\nHalf-mass radius: " << Rhm
              << "\nCircular Velocity: " << vc << "\nMean interparticle separation: " << mips << "\ntcross: " << tcross
              << "\ntrelax: " << trelax << std::endl;

    // Export stats
    std::ofstream stats("../results/test_stats.txt");
    stats << particles[N - 1].radius() << std::endl
          << quarter_mass_r << std::endl
          << Rhm << std::endl
          << vc << std::endl
          << tcross << std::endl
          << trelax << std::endl;

    std::ofstream results("../results/radii.txt");
    for (int i = 0; i < N; ++i)
        results << particles[i].radius() << " ";
    results << std::endl;
}

// calculate the density function by inferring it from the particle distribution as well as the analytical Hernquist density function
void calc_distr(std::vector<Particle> &particles, double M, double quarter_mass_r)
{
    int N = std::min(particles.size(), size_t{NUM_PARTICLES});
    // Count occurences in shells with equal number of particles
    double r = 0.;
    double mass = 0.;
    double volume;
    double single_mass = M / double{N};
    double hernquist_density;
    double r_old = 0;
    std::vector<double> densities, shells, hernquist, stdev;

    int bin = 1;
    int particles_per_shell = 1000;

    // we want the volume of the shell, not the sphere, and hernquist measures
    // somewhere in the middle of the shell, so average radius
    for (int i = 0; i < N; ++i)
    {
        mass += particles[i].mass();
        if ((i >= bin * particles_per_shell) || (i == N - 1))
        {
            r = particles[i].radius();
            volume = (4. / 3. * M_PI * (r * r * r - r_old * r_old * r_old));
            hernquist_density = M / (2. * M_PI) * quarter_mass_r / ((r + r_old) / 2 * std::pow((r + r_old) / 2 + quarter_mass_r, 3));
            densities.push_back(mass / volume);
            shells.push_back((r + r_old) / 2);
            hernquist.push_back(hernquist_density);
            stdev.push_back(std::sqrt(hernquist_density * volume / single_mass) * single_mass / volume);
            mass = 0;
            r_old = r;
            ++bin;
        }
    }

    // std::cout << shells.size() << " " << densities.size() << " " << hernquist.size() << " " << stdev.size() << std::endl;

    // Export results for plotting
    std::ofstream results("../results/test_distr.txt");
    int end = densities.size();
    for (int i = 0; i < end; ++i)
        results << densities[i] << " ";
    results << std::endl;
    for (int i = 0; i < end; ++i)
        results << shells[i] << " ";
    results << std::endl;
    for (int i = 0; i < end; ++i)
        results << hernquist[i] << " ";
    results << std::endl;
    for (int i = 0; i < end; ++i)
        results << hernquist[i] + stdev[i] << " ";
    results << std::endl;
    for (int i = 0; i < end; ++i)
        results << hernquist[i] - stdev[i] << " ";
}

// calculate the forces acting upon one particle
Vector p_force(std::vector<Particle> &particles, Vector position, double eps, const int length, int id)
{
    Vector force, dif;
    double denom;
    force = Vector::Zero();

    for (int i = 0; i < length; ++i)
    {
        dif = position - particles[i].position();
        denom = std::pow(dif.squaredNorm() + eps * eps, 1.5);
        if (i != id)
            force += -G * particles[i].mass() / denom * dif;
    }
    return force;
}

// calculate the forces on each particle using a brute-force n-body approach
std::vector<Vector> direct_force_calculation(std::vector<Particle> &particles, double eps)
{

    int length = std::min(particles.size(), size_t{NUM_PARTICLES});

    std::vector<Vector> forces(length, Vector::Zero());
    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; ++i)
    {
        if (i == 0)
            std::cout << "Num threads: " << omp_get_num_threads() << std::endl;
        forces[i] = particles[i].mass() * p_force(particles, particles[i].position(), eps, length, i);
    }

    /*
    unsigned int length = particles.size();
    std::vector<Vector> forces(length, Vector::Zero());
    Vector force, dif;
    double denom;
    auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < length; ++i){
        //if (i % 1000 == 0) std::cout << "Iteration " << i << " reached\n";
        mass = particles[i].mass();
        for (unsigned int j = i+1; j < length; ++j){
            dif = particles[i].position() - particles[j].position();
            denom = std::pow(dif.squaredNorm() + eps*eps, 1.5);
            force = -G * mass * partiparticles[j].mass() / denom * dif;
            forces[j] -= force;
            forces[i] += force;
        }
    }
    */
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout << "Time for force calculation: " << time.count() << "s" << std::endl;

    return forces;
}

// calculate the analytical forces expected by Newton's second theorem for spherical potentials
void newton_forces(std::vector<Particle> particles)
{
    int length = std::min(particles.size(), size_t{NUM_PARTICLES});

    std::vector<Vector> forces;
    std::vector<double> shells;
    Vector r;
    double rnorm, r3;
    double mass = 0;

    for (int i = 0; i < length; ++i)
    {
        mass += particles[i].mass();
        if (!(i % 1))
        {
            r = particles[i].position();
            rnorm = r.norm();
            r3 = rnorm * rnorm * rnorm;
            forces.push_back(-G * particles[i].mass() * mass / r3 * r);
            shells.push_back(rnorm);
        }
    }

    std::ofstream results("../results/newton_forces.txt");
    double its = forces.size();
    for (int i = 0; i < its; ++i)
        results << shells[i] << " ";
    results << std::endl;
    for (int i = 0; i < its; ++i)
        results << forces[i].norm() << " ";
    results << std::endl;
}

// function to print the forces to output files
void print_state(std::vector<Particle> &particles, std::vector<Vector> &forces, std::ofstream &results)
{
    int length = std::min(particles.size(), size_t{NUM_PARTICLES});

    for (int i = 0; i < length; ++i)
        results << std::setw(15) << forces[i].norm() * units[3] << " ";
    results << std::endl;
    for (int i = 0; i < length; ++i)
        results << std::setw(15) << particles[i].x() * units[0] << " ";
    results << std::endl;
    for (int i = 0; i < length; ++i)
        results << std::setw(15) << particles[i].y() * units[0] << " ";
    results << std::endl;
    for (int i = 0; i < length; ++i)
        results << std::setw(15) << particles[i].z() * units[0] << " ";
    results << std::endl;
}

int main()
{
    std::ifstream data("../data/data.txt");
    std::vector<Particle> particles;
    double M = 0;
    double quarter_mass_r;
    double tcross;
    double dt = 1e-6; // FIND CORRECT VALUE

#pragma omp master
    {
        std::cout << "OMP max threads: " << omp_get_max_threads() << std::endl;
    }

    std::cout << "\nDetermining the units\n";
    det_units();

    // Read the data, create vectors of the radii and the masses and calculate total mass
    std::cout << "\nPreparing the data\n";
    prep_data(data, particles, M, quarter_mass_r, tcross);

    // Calculate the distribution profile of the data and compare it to theory
    std::cout << "\nCalculating the distribution\n";
    calc_distr(particles, M, quarter_mass_r);

    newton_forces(particles);

    // Calculate the direct forces using different eps and integrate them
    std::cout << "Calculating the forces\n";
    // int length = particles.size();
    int length = std::min(particles.size(), size_t{NUM_PARTICLES});

    std::vector<Vector> forces(length, Vector::Zero());
    std::vector<Vector> approx_forces(length, Vector::Zero());
    Vector vnew, xnew;
    std::vector<double> eps = {0.1}; //{0.0, 0.01, 0.1, 1.0};
    std::vector<std::string> folders = {"../results/brute_test_run/"};
    //{"../results/integration_0.01_1e-6/"};
    //{"../results/integration_0.0/", "../results/integration_0.01/", results/integration_0.1/", "../results/integration_1.0/"};
    for (int e = 0; e < eps.size(); ++e)
    {
        std::cout << "Starting simulation for eps = " << eps[e] << std::endl;
        std::ofstream out(folders[e] + "t_0.000000.txt");

        print_state(particles, forces, out);
        // for (double t = 0; t < 10 * tcross; t += dt){
        for (double t = 0; t < 1; t = 1.1)
        {
            forces = direct_force_calculation(particles, eps[e]);
            for (int i = 0; i < length; ++i)
            { // simple symplectic Euler, for now
                Particle *p = &particles[i];
                vnew = p->velocity() + dt * forces[i] / particles[i].mass();
                xnew = p->position() + dt * vnew;
                p->set_position(xnew);
                p->set_velocity(vnew);
            }
            std::cout << "Timestep " << t << " complete\n";
            std::ofstream out(folders[e] + "t_" + std::to_string(t + dt) + ".txt");
            print_state(particles, forces, out);
        }
    }

    // Build the oct-tree
    std::cout << "\nBuilding the oct-tree\n";
    OctTree tree(length);
    tree.BuildTree(particles);

    // Calculate forces on the oct-tree with opening angle theta_c
    std::cout << "Calculating the forces on the tree\n";
    std::vector<double> thetas = {0.02}; //{0., 0.001, 0.01, 0.1, 1.};
    for (int i = 0; i < thetas.size(); ++i)
    {
        std::string fileName = "../results/tree_test_run.txt"; // "../results/tree_forces/parallel" + std::to_string(thetas[i]) + ".txt";
        std::ofstream out(fileName);
        approx_forces = tree.CalculateTreeForces(particles, thetas[i], G);
        print_state(particles, approx_forces, out);
    }

    return 0;
}
