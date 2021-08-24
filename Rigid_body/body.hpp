//
// Created by jacob on 31/7/21.
//

#ifndef CODE_BODY_HPP
#define CODE_BODY_HPP

#include <vector>
//#include <memory>
#include "../MyMath/vec3.hpp"


#include "../MyMath/ops.hpp"
#include <numeric>
#include <functional>
#include <fstream>
#include "mesh.hpp"

struct body {
    mesh model; //positions and velocities global

    vec3 pos_cm;  //the position of the center of mass
    vec3 vel_cm{};  //the velocity of the center of mass

    const double M;   //the mass of the system

    vec3 w_cm{};    //the angular velocity of the system about the center of mass

    body() = delete;
    body(const std::vector<vec3> &pos_, const std::vector<double>& mass_) : model(pos_, mass_), M(std::accumulate(mass_.begin(), mass_.end(), 0.0)),
                                         pos_cm{ std::inner_product(mass_.begin(), mass_.end(), pos_.begin(), vec3{}) / std::accumulate(mass_.begin(), mass_.end(), 0.0)}{
#ifndef NDEBUG
        if (pos_.size() != mass_.size()) {
            std::cerr << "mass and pos must be the same size\n";
        }
#endif
    }

    explicit body(const mesh &model_) : model(model_), M(std::accumulate(model_.mass.begin(), model_.mass.end(), 0.0)),
            pos_cm{ std::inner_product(model_.mass.begin(), model_.mass.end(), model_.vertices.begin(), vec3{}) / std::accumulate(model_.mass.begin(), model_.mass.end(), 0.0)} {}




    //finds the positions of all the particles undergoing forces <forces> after time dt
    void update_pos(const std::vector<vec3> &forces, const double dt) {
#ifndef DEBUG
        if (forces.size() != model.vertices.size()) {
            std::cerr << "Need forces for every vertex\n";
        }
#endif

        //finding the linear acceleration
        const vec3 a =  std::accumulate(forces.begin(), forces.end(), vec3{}) / M;

        //finding position and velocity of CoM
        vel_cm += a*dt;
        const auto pos_cm_old = pos_cm;
        pos_cm += vel_cm*dt;

        vec3 test{};
        for (unsigned i = 0; i < forces.size(); i++) {
            test += cross(model.vertices[i] - pos_cm, forces[i]);
            //std::cout << test << "\t" << (*pos)[i] << "\t" << pos_cm << "\t" << (*pos)[i] - pos_cm << "\t" << forces[i] << "\n";
        }

        //finding the torque
        const vec3 t = std::inner_product(model.vertices.begin(), model.vertices.end(), forces.begin(), vec3{}, std::plus<>(),
                                          [&](const vec3 &r, const vec3 &F){return cross( (r-pos_cm), F);});

        //moment of inertia
        const c_line r_axis(pos_cm, t);
        const double I = std::inner_product(model.mass.begin(), model.mass.end(), model.vertices.begin(), 0.0, std::plus<>(),
                                            [&](const double &m, const vec3 &r){return m * r_axis.distance(r)*r_axis.distance(r);});

        //angular acceleration
        const vec3 alpha = t/I;

        //finding angle and angular-velocity pseudo-vector
        w_cm += alpha*dt;
        const vec3 rot_angle_vec = w_cm*dt;

        //updating all positions
        std::transform(model.vertices.begin(), model.vertices.end(), model.vertices.begin(),
                       [&](const vec3& x)
                       {return rotate(pos_cm_old, rot_angle_vec, x, rot_angle_vec.length()) + vel_cm*dt;}); //rotating about the old center of mass, then moving forward

        //updating all velocities
        std::transform(model.velocities.begin(), model.velocities.end(), model.velocities.begin(),
                       [&](const vec3& x)
                       {return vel_cm + cross( (x-pos_cm),w_cm);}); //velocity of cm plus rotational velocity (v=r x w)
                                                //https://en.wikipedia.org/wiki/Angular_velocity
    }



    void debug_write_pos(const char* file_loc) {
        std::ofstream output(file_loc);
        for (const auto& p : model.vertices) {
            output << p << "\n";
        }
        output.close();
    }
};







#endif //CODE_BODY_HPP
