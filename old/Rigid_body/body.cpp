//
// Created by jacob on 31/7/21.
//

#include "body.hpp"
#include "../MyMath/ops.hpp"
#include <numeric>
#include <functional>
#include <fstream>

body::body(std::shared_ptr<std::vector<vec3>> &pos_, const std::vector<double>& mass_) : pos(pos_), mass(mass_),
    M(std::accumulate(mass_.begin(), mass_.end(), 0.0)),
    pos_cm{ std::inner_product(mass_.begin(), mass_.end(), pos_->begin(), vec3{}) / std::accumulate(mass_.begin(), mass_.end(), 0.0)}{
#ifndef NDEBUG
    if (pos_->size() != mass.size()) {
        std::cerr << "mass and pos struct must be the same size\n";
    }
#endif
}

void body::update_pos(const std::vector<vec3> &forces, const double dt) {
    //finding the linear acceleration
    const vec3 a =  std::accumulate(forces.begin(), forces.end(), vec3{}) / M;

    //finding position and velocity of CoM
    vel_cm += a*dt;
    const auto pos_cm_old = pos_cm;
    pos_cm += vel_cm*dt;
    //std::cout << vel_cm << "\n";
    //std::cout << pos_cm << "\n";

    //for (const auto& f : forces) {
    //    std::cout << f << "\n";
    //}

    vec3 test{};
    for (unsigned i = 0; i < forces.size(); i++) {
        test += cross((*pos)[i] - pos_cm, forces[i]);
        //std::cout << test << "\t" << (*pos)[i] << "\t" << pos_cm << "\t" << (*pos)[i] - pos_cm << "\t" << forces[i] << "\n";
    }

    //finding the torque
    const vec3 t = std::inner_product(pos->begin(), pos->end(), forces.begin(), vec3{}, std::plus<>(),
            [&cm = pos_cm](const vec3 &r, const vec3 &F){return cross( (r-cm), F);});

    //moment of inertia
    const c_line r_axis(pos_cm, t);
    const double I = std::inner_product(mass.begin(), mass.end(), pos->begin(), 0.0, std::plus<>(),
                              [&l = r_axis](const double &m, const vec3 &r){return m * l.distance(r)*l.distance(r);});

    //angular acceleration
    const vec3 alpha = t/I;

    //finding angle and angular-velocity pseudo-vector
    w_cm += alpha*dt;
    const vec3 rot_angle_vec = w_cm*dt;
    //std::cout << I << "\n";
    //std::cout << "\n" << alpha << "\t" << a << "\t" << t << "\t" << rot_angle_vec << "\n";

    //std::cout << (*pos)[12];
    //updating all positions
    std::transform(pos->begin(), pos->end(), pos->begin(),
                   [d = vel_cm*dt, theta = rot_angle_vec, cm = pos_cm_old](const vec3& x)
                   {return rotate(cm, theta, x, theta.length()) + d;});
    //std::cout << vel_cm*dt << "\t" << rotate(pos_cm_old, rot_angle_vec, (*pos)[12], rot_angle_vec.length()) << "\n";
    //std::cout << "\t" << (*pos)[12] << "\n";
}


void body::debug_write_pos(const char* file_loc) {
    std::ofstream output(file_loc);
    for (const auto& p : *pos) {
        output << p << "\n";
    }
}