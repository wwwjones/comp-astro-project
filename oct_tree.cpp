#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Core>
#include "particle.hpp"
#include "cell.hpp"
#include "oct_tree.hpp"

void OctTree::BuildTree(std::vector<Particle> &particles)
{
    // Place particles one at a time into the tree
    int length = std::min(particles.size(), max_particles_);

    double sidelength = 1.;
    double max_radius = particles.back().radius();
    while (sidelength < 2. * max_radius)
    {
        sidelength *= 2.;
    }

    tree_ = new Cell(0, particles[0].position(), particles[0].mass());
    ++num_particles_;
    ++num_cells_;

    // Build the tree
    for (int i = 1; i < length; ++i)
    {
        level_ = 0;
        Place(i, particles, Vector::Zero(), sidelength, tree_);
        if (level_ > max_level_)
            max_level_ = level_;
    }

    // Calculate the quadrupoles
    CalculateProperties(tree_);

    std::cout << "Tree fully built" << std::endl;
    std::cout << "Particles placed: " << num_particles_ << std::endl;
    std::cout << "Max sidelength: " << sidelength << std::endl;
    std::cout << "Center of mass: \n"
              << tree_->com() << std::endl;
    std::cout << "Number of levels: " << max_level_ << std::endl;
    std::cout << "Number of cubes: " << num_cells_ << std::endl;
}

// Place a particle into the tree so that each particle is the sole inhabitant of a cell
void OctTree::Place(int i, std::vector<Particle> &particles, Vector center, double sidelength, Cell *current)
{
    ++level_;
    Particle p = particles[i];
    Vector position = p.position();
    double mass = p.mass();
    int id = current->id();
    std::pair<int, Vector> idx_cntr; // for storing the index and the new center
    idx_cntr = PickNextCell(position, center, sidelength / 4.);
    int index = idx_cntr.first;
    Vector next_center = idx_cntr.second;
    // Update this cubes info as we are passing through
    current->update_cell(position, mass); // update the mass and CoM

    if (id >= 0)
    { // If the current cell has an occupant
        // Put our current particle in the right cell
        Cell *next_cube = new Cell(i, position, mass);
        ++num_cells_;
        current->set_next(next_cube, index);
        // We have to re-place the inhabitant of this cell
        current->rmv_from_cell(particles[id]); // reverse the CoM and mass because it will be r>
        current->set_id(-1);
        --level_;
        Place(id, particles, center, sidelength, current);
    }
    else
    {
        // If the next cell exists, continue tree traversal, else create a new cell and place the particle
        Cell *next = current->next(index);
        if (next)
        {
            current = next;
            Place(i, particles, next_center, sidelength / 2., current);
        }
        else
        {
            Cell *next_cube = new Cell(i, position, mass);
            ++num_particles_;
            ++num_cells_;
            current->set_next(next_cube, index);
            ++level_;
        }
    }
}

// Returns the index for the correct next cell pointer and updates the center
std::pair<int, Vector> OctTree::PickNextCell(Vector position, Vector center, double sidelength)
{
    // returns the id of the next cell
    // maybe a shorter implementation? center += (-1 ^ (position < center)) * sidelength
    std::pair<int, Vector> cell_and_center;
    if (position.x() < center.x())
    {
        if (position.y() < center.y())
        {
            if (position.z() < center.z())
            {
                center.x() -= sidelength;
                center.y() -= sidelength;
                center.z() -= sidelength;
                cell_and_center.first = 0;
                cell_and_center.second = center;
                return cell_and_center; // dlf
            }
            else
            {
                center.x() -= sidelength;
                center.y() -= sidelength;
                center.z() += sidelength;
                cell_and_center.first = 1;
                cell_and_center.second = center;
                return cell_and_center; // ulf
            }
        }
        if (position.z() < center.z())
        {
            center.x() -= sidelength;
            center.y() += sidelength;
            center.z() -= sidelength;
            cell_and_center.first = 2;
            cell_and_center.second = center;
            return cell_and_center; // dlb
        }
        else
        {
            center.x() -= sidelength;
            center.y() += sidelength;
            center.z() += sidelength;
            cell_and_center.first = 3;
            cell_and_center.second = center;
            return cell_and_center; // ulb
        }
    }
    if (position.y() < center.y())
    {
        if (position.z() < center.z())
        {
            center.x() += sidelength;
            center.y() -= sidelength;
            center.z() -= sidelength;
            cell_and_center.first = 4;
            cell_and_center.second = center;
            return cell_and_center; // drf
        }
        else
        {
            center.x() += sidelength;
            center.y() -= sidelength;
            center.z() += sidelength;
            cell_and_center.first = 5;
            cell_and_center.second = center;
            return cell_and_center; // urf
        }
    }
    if (position.z() < center.z())
    {
        center.x() += sidelength;
        center.y() += sidelength;
        center.z() -= sidelength;
        cell_and_center.first = 6;
        cell_and_center.second = center;
        return cell_and_center; // drb
    }
    else
    {
        center.x() += sidelength;
        center.y() += sidelength;
        center.z() += sidelength;
        cell_and_center.first = 7;
        cell_and_center.second = center;
        return cell_and_center; // urb
    }
}

// Recursively calculate the quadrupole and the radius of each cell
void OctTree::CalculateProperties(Cell *current)
{
    int n = 0;
    // calculate the quadrupole and count how many particles are in the cell
    Matrix Q = CalculateQuadrupole(current->com(), n, current);
    current->set_quadrupole(Q);
    // calculate the average distance from the com
    double r = CalculateRadius(current->com(), n, current);
    current->set_radius(r);

    for (int i = 0; i < 8; ++i)
    {
        Cell *next = current->next(i);
        if (next)
            CalculateProperties(next);
    }
}

// Recursively calculate a cell's quadrupole with respect to a given center of mass
Matrix OctTree::CalculateQuadrupole(Vector com, int &n, Cell *current)
{
    Matrix Q = Matrix::Zero();
    if (current->id() >= 0)
    {
        ++n;
        Vector dist = com - current->com();
        for (int j = 0; j < 3; ++j)
        {
            for (int i = j; i < 3; ++i)
            {
                Q(i, j) = current->mass() * (3 * dist[i] * dist[j] - (i == j) * dist.squaredNorm());
            }
        }
        Q(0, 1) = Q(1, 0);
        Q(0, 2) = Q(2, 0);
        Q(1, 2) = Q(2, 1);
    }

    else
    {
        for (int i = 0; i < 8; ++i)
        {
            Cell *next = current->next(i);
            if (next)
            {
                Q += CalculateQuadrupole(com, n, current->next(i));
            }
        }
    }
    return Q;
}

// Recursively calculate a cell's radius with respect to a given center of mass
double OctTree::CalculateRadius(Vector com, int n, Cell *current)
{
    double r = 0;
    if (current->id() >= 0)
    {
        r = (current->com() - com).norm();
        r /= n;
    }

    else
    {
        for (int i = 0; i < 8; ++i)
        {
            Cell *next = current->next(i);
            if (next)
            {
                r += CalculateRadius(com, n, current->next(i));
            }
        }
    }
    return r;
}

std::vector<Vector> OctTree::CalculateTreeForces(std::vector<Particle> &particles, double theta_c, double G)
{
    int length = std::min(particles.size(), max_particles_);

    std::vector<Vector> approx_forces(length, Vector::Zero());

    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; ++i)
    {
        approx_forces[i] = particles[i].mass() * CalculateForces(i, particles[i].position(), theta_c, G, tree_);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout << "Time for approximate force calculation: " << time.count() << "s" << std::endl;

    return approx_forces;
}

Vector OctTree::CalculateForces(int id, Vector position, double theta_c, double G, Cell *current)
{
    Vector force = Vector::Zero();
    int cur_id = current->id();
    if (cur_id == id)
        return force; // Don't calculate self force, divides by 0
    Vector y = position - current->com();
    double y_magn = y.norm();
    double theta = current->radius() / y_magn;
    theta = 2.0 * std::asin(theta);
    if (theta > 1)
        theta = M_PI;
    if ((theta <= theta_c))
    {
        Matrix Q = current->Q();
        double M = current->mass();
        double temp, yi;
        double y2 = y_magn * y_magn;
        double y3 = y2 * y_magn;
        double y7 = y3 * y2 * y2;
        double yQy = y.transpose() * Q * y;
        for (int i = 0; i < 3; ++i)
        {
            yi = y(i);
            temp = 0;
            for (int j = 0; j < 3; ++j)
            {
                temp += Q(i, j) * y(j);
            }
            force(i) = -G * ((-M * yi / y3) + ((temp * y2 - yQy * 2.5 * yi) / y7));
        }
    }

    else
    {
        for (int i = 0; i < 8; ++i)
        {
            Cell *next = current->next(i);
            if (next)
                force += CalculateForces(i, position, theta_c, G, next);
        }
    }
    return force;
}