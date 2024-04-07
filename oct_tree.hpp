#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Core>
#include "particle.hpp"
#include "cell.hpp"

#ifndef SIMULATION_H
#define SIMULATION_H

class OctTree
{
public:
    // Default constructor
    OctTree() : tree_(new Cell), max_particles_(UINT_MAX), num_particles_(0), num_cells_(0), max_level_(0), level_(0), side_length_(0), com_(Eigen::Vector3d::Zero()) {}

    // Constructor with a given maximum number of particles to process
    OctTree(int max_particles) : tree_(new Cell), max_particles_(max_particles), num_particles_(0), num_cells_(0), max_level_(0), level_(0), side_length_(0), com_(Eigen::Vector3d::Zero()) {}

    // Build the oct-tree by placing the particles into the tree sequentially then calculate the tree's properties recursively
    void BuildTree(std::vector<Particle> &particles);

    // Calculate the forces acting upon all the particles in the tree using a given opening angle
    std::vector<Vector> CalculateTreeForces(std::vector<Particle> &particles, double theta_c, double G);

private:
    Cell *tree_;           // the head of the oct-tree
    size_t max_particles_; // the maximum number of particles we want to process
    int num_particles_;    // the number of particles placed in the tree
    int num_cells_;        // the number of cells in the tree
    int max_level_;        // the depth of the tree
    int level_;            // the current depth of the tree
    int side_length_;      // the side length of the cubic bounding box of all particles in the system
    Eigen::Vector3d com_;  // the center of mass of the tree

    // Place a particle into the tree so that each particle is the sole inhabitant of a cell
    void Place(int i, std::vector<Particle> &particles, Vector center, double sidelength, Cell *current);

    // Returns the index for the correct next cell pointer and updates the center
    std::pair<int, Vector> PickNextCell(Vector position, Vector center, double sidelength);

    // Recursively calculate the quadrupole and the radius of each cell
    void CalculateProperties(Cell *current);

    // Recursively calculate a cell's quadrupole with respect to a given center of mass
    Matrix CalculateQuadrupole(Vector com, int &n, Cell *current);

    // Recursively calculate a cell's radius with respect to a given center of mass
    double CalculateRadius(Vector com, int n, Cell *current);

    // Calculate the forces acting upon one particle
    Vector CalculateForces(int id, Vector position, double theta_c, double G, Cell *current);
};

#endif