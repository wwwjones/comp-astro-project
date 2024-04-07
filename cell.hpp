#include <vector>
#include <Eigen/Core>
#include "particle.hpp"

#ifndef MULTIPOLE_H
#define MULTIPOLE_H

// int total_p = 0; // counts particles placed
// int total_c = 0; // counts number of cubes in tree

using Matrix = Eigen::Matrix3d;
using Vector = Eigen::Vector3d;

class Cell
{
public:
    Cell() : id_(-1), com_(Vector::Zero()), radius_(0.),
             mass_(0.), quadrupole_(Matrix::Zero())
    {
        //++total_p;
        //++total_c;
        for (int i = 0; i < 8; ++i)
            next_.push_back(nullptr);
    }
    Cell(int id, Vector position, double mass) : id_(id), com_(position), radius_(position.norm()),
                                                 mass_(mass), quadrupole_(Matrix::Zero())
    {
        //++total_p;
        //++total_c;
        for (int i = 0; i < 8; ++i)
            next_.push_back(nullptr);
    }
    // setters
    void set_id(int id) { id_ = id; }
    void set_radius(double r) { radius_ = r; }
    void set_quadrupole(Matrix Q) { quadrupole_ = Q; }
    void set_next(Cell *next, int i) { next_[i] = next; }
    // getters
    int id() { return id_; }
    Vector com() { return com_; }
    double radius() { return radius_; }
    double mass() { return mass_; }
    Matrix Q() { return quadrupole_; }
    Cell *next(int i) { return next_[i]; }
    /// other member fucntions
    void update_cell(Vector position, double m)
    {
        com_ *= mass_;
        com_ += m * position;
        mass_ += m;
        com_ /= mass_;
    }
    void rmv_from_cell(Particle p)
    {
        com_ *= mass_;
        com_ -= p.mass() * p.position();
        mass_ -= p.mass();
        com_ /= mass_;
        //--total_p;
    }
    ~Cell()
    {
        id_ = radius_ = mass_ = 0;
        com_ = Vector::Zero();
        quadrupole_ = Matrix::Zero();
        // also call delete on all the ptrs and set them to 0
        for (int i = 0; i < 8; ++i)
        {
            if (next_[i])
                delete next_[i];
            next_[i] = nullptr;
        }
    }

private:
    int id_;                   // id of the particle that is in it, -1 if empty
    Vector com_;               // position of the sector's center of mass
    double radius_;            //
    double mass_;              // total mass of the sector
    Matrix quadrupole_;        // quadrupole moment of the cell's group of particles
    std::vector<Cell *> next_; // vector of this cells children
    //[dlf, ulf, dlb, uln, drf, urf, drb, urb]
    //[ 0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ]
    // up/down, left/right, front/back
};

#endif
