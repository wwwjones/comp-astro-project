#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "particle.h"

#ifndef MULTIPOLE_H
#define MULTIPOLE_H

unsigned int total_p = 0; //counts particles placed
unsigned int total_c = 0; //counts number of cubes in tree

using Matrix = Eigen::Matrix3d;
using Vector = Eigen::Vector3d;

class cube {
    public:
        cube() : id_(-1), com_(Vector::Zero()), radius_(0.),
            mass_(0.), quadrupole_(Matrix::Zero()) {
            ++total_p;
            ++total_c;
            for (int i = 0; i < 8; ++i)
                next_.push_back(nullptr);
            }
        cube(particle p, int id, Vector center, double sidelength) :
            id_(id), com_(p.position()), radius_(p.position().norm()),
            mass_(p.mass()), quadrupole_(Matrix::Zero()) {
            ++total_p;
            ++total_c;
            //std::cout << id << " was placed around\n" << center_ << std::endl;
            for (int i = 0; i < 8; ++i)
                next_.push_back(nullptr);
            }
        //setters
        void set_id(int id) {id_ = id;}
        void set_radius(double r) {radius_ = r;}
        void set_quadrupole(Matrix Q) {quadrupole_ = Q;}
        void set_next(cube* next, int i) {next_[i] = next;}
        //getters
        int id() {return id_;}
        Vector com() {return com_;}
        double radius() {return radius_;}
        Matrix Q() {return quadrupole_;}
        cube* next(int i) {return next_[i];}
        ///other member fucntions
        void update_cube(particle p){
            com_ *= mass_;
            com_ += p.mass() * p.position();
            mass_ += p.mass();
            com_ /= mass_;
        }
        void rmv_from_cube(particle p){
            com_ *= mass_;
            com_ -= p.mass() * p.position();
            mass_ -= p.mass();
            com_ /= mass_;
            --total_p;
        }
        //void update_radius(){};
        //void update_quadrupole(){};
        ~cube() {
            id_ = radius_ = mass_ = 0;
            com_ = Vector::Zero();
            quadrupole_ = Matrix::Zero();
            //also call delete on all the ptrs and set them to 0
            for (int i = 0; i < 8; ++i) {
                if (next_[i]) delete next_[i];
                next_[i] = nullptr;
            }
        }
    private:
        int id_; //id of the particle that is in it, -1 if empty
        Vector com_;
        double radius_;
        double mass_;
        Matrix quadrupole_;
        std::vector<cube*> next_;
        //[dlf, ulf, dlb, uln, drf, urf, drb, urb]
        //[ 0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ]
        // up/down, left/right, front/back
};

#endif
