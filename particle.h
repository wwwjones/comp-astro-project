#include <iostream>
#include <vector>
#include <Eigen/Core>

#ifndef PARTICLE_H
#define PARTICLE_H

class particle {
    public:
        particle() : position_(Eigen::Vector3d::Zero()), mass_(0.) {}
        //setters
        void set_x(double x) {position_.x() = x;}
        void set_y(double y) {position_.y() = y;}
        void set_z(double z) {position_.z() = z;}
        void set_mass(double m) {mass_ = m;}
        //getters
        double x() {return position_.x();}
        double y() {return position_.y();}
        double z() {return position_.z();}
        Eigen::Vector3d position() {return position_;}
        double mass() {return mass_;}
    private:
        Eigen::Vector3d position_;
        double mass_;
};

#endif
