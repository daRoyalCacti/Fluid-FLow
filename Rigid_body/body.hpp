//
// Created by jacob on 31/7/21.
//

#ifndef CODE_BODY_HPP
#define CODE_BODY_HPP

#include <vector>
#include <memory>
#include "../MyMath/vec3.hpp"


#include "../MyMath/ops.hpp"
#include <numeric>
#include <functional>
#include <fstream>

struct body {
    std::shared_ptr<std::vector<vec3>> pos;    //the position of every particle
    const std::vector<double> mass;    //the mass of every particle

    vec3 pos_cm;  //the position of the center of mass
    vec3 vel_cm{};  //the velocity of the center of mass

    const double M;   //the mass of the system

    vec3 w_cm{};    //the angular velocity of the system about the center of mass

    body() = delete;
    body(std::shared_ptr<std::vector<vec3>> &pos_, const std::vector<double>& mass_) : pos(pos_), mass(mass_),
                                                                                       M(std::accumulate(mass_.begin(), mass_.end(), 0.0)),
                                                                                       pos_cm{ std::inner_product(mass_.begin(), mass_.end(), pos_->begin(), vec3{}) / std::accumulate(mass_.begin(), mass_.end(), 0.0)}{
#ifndef NDEBUG
        if (pos_->size() != mass.size()) {
            std::cerr << "mass and pos struct must be the same size\n";
        }
#endif
    }

    //finds the positions of all the particles undergoing forces <forces> after time dt
    void update_pos(const std::vector<vec3> &forces, const double dt) {
        //finding the linear acceleration
        const vec3 a =  std::accumulate(forces.begin(), forces.end(), vec3{}) / M;

        //finding position and velocity of CoM
        vel_cm += a*dt;
        const auto pos_cm_old = pos_cm;
        pos_cm += vel_cm*dt;

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

        //updating all positions
        std::transform(pos->begin(), pos->end(), pos->begin(),
                       [d = vel_cm*dt, theta = rot_angle_vec, cm = pos_cm_old](const vec3& x)
                       {return rotate(cm, theta, x, theta.length()) + d;});
    }



    void debug_write_pos(const char* file_loc) {
        std::ofstream output(file_loc);
        for (const auto& p : *pos) {
            output << p << "\n";
        }
        output.close();
    }
};







#endif //CODE_BODY_HPP
