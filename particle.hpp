#include <iostream>
#include <vector>
#include <Eigen/Core>

#ifndef PARTICLE_H
#define PARTICLE_H

class Particle
{
public:
    Particle() : position_(Eigen::Vector3d::Zero()), velocity_(Eigen::Vector3d::Zero()), radius_(0.), mass_(0.) {}
    // setters
    void set_x(double x) { position_.x() = x; }
    void set_y(double y) { position_.y() = y; }
    void set_z(double z) { position_.z() = z; }
    void set_position(Eigen::Vector3d x) { position_ = x; }
    void set_vx(double vx) { velocity_.x() = vx; }
    void set_vy(double vy) { velocity_.y() = vy; }
    void set_vz(double vz) { velocity_.z() = vz; }
    void set_velocity(Eigen::Vector3d v) { velocity_ = v; }
    void set_radius(double r) { radius_ = r; }
    void set_mass(double m) { mass_ = m; }
    // getters
    double x() { return position_.x(); }
    double y() { return position_.y(); }
    double z() { return position_.z(); }
    Eigen::Vector3d position() { return position_; }
    double vx() { return velocity_.x(); }
    double vy() { return velocity_.y(); }
    double vz() { return velocity_.z(); }
    Eigen::Vector3d velocity() { return velocity_; }
    double radius() { return radius_; }
    double mass() { return mass_; }

private:
    Eigen::Vector3d position_;
    Eigen::Vector3d velocity_;
    double radius_;
    double mass_;
};

#endif
