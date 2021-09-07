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
    // - points are the positions where the foces are applied
    void update_pos(const std::vector<vec3> &forces, const std::vector<vec3> & points, const double dt) {
#ifndef DEBUG
        if (forces.size() != points.size()) {
            std::cerr << "Need forces for every vertex\n";
        }
#endif

        //finding the linear acceleration
        const vec3 a =  std::accumulate(forces.begin(), forces.end(), vec3{}) / M;

        //finding position and velocity of CoM
        vel_cm += a*dt;
        const auto pos_cm_old = pos_cm;
        pos_cm += vel_cm*dt;

        //std::cerr << vel_cm << "\t" << pos_cm << "\n";


        //finding the torque
        const vec3 t = std::inner_product(points.begin(), points.end(), forces.begin(), vec3{}, std::plus<>(),
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

#ifndef NDEBUG
        {
            bool err = false;
            if (!std::isfinite(vel_cm.x()) || !std::isfinite(vel_cm.y()) || !std::isfinite(vel_cm.z())) {
                err = true;
                std::cerr << "rigid body got infinite velocity\n";
            }
            if (!std::isfinite(pos_cm.x()) || !std::isfinite(pos_cm.y()) || !std::isfinite(pos_cm.z())) {
                err = true;
                std::cerr << "rigid body got infinite position\n";
            }
            if (!std::isfinite(t.x()) || !std::isfinite(t.y()) || !std::isfinite(t.z())) {
                err = true;
                std::cerr << "rigid body got infinite torque\n";
            }
            if (!std::isfinite(I)) {
                std::cerr << "rigid body got infinite moment of inertia\n";
                err = true;
            }
            if (!std::isfinite(alpha.x()) || !std::isfinite(alpha.y()) || !std::isfinite(alpha.z())) {
                err = true;
                std::cerr << "rigid body got infinite angular acceleration\n";
            }
            if (!std::isfinite(w_cm.x()) || !std::isfinite(w_cm.y()) || !std::isfinite(w_cm.z())) {
                err = true;
                std::cerr << "rigid body got infinite angular velocity\n";
            }
            if (!std::isfinite(rot_angle_vec.x()) || !std::isfinite(rot_angle_vec.y()) || !std::isfinite(rot_angle_vec.z())) {
                err = true;
                std::cerr << "rigid body got infinite rotation angle\n";
            }

            if (err) {
                std::cerr << "\tpos : " <<pos_cm << "\n";
                std::cerr << "\told pos : " << pos_cm_old << "\n";
                std::cerr << "\tvelocity : " << vel_cm << "\n";
                std::cerr << "\ttorque : " << t << "\n";
                std::cerr << "\tmoment of inertia : " << I << "\n";
                std::cerr << "\tangular accelerator : " << alpha << "\n";
                std::cerr << "\tangular velocity : " << w_cm << "\n";
                std::cerr << "\trotation vector : " << rot_angle_vec << "\n";
            }
        }
#endif





        //updating all positions
        std::transform(model.vertices.begin(), model.vertices.end(), model.vertices.begin(),
                       [&](const vec3& x)
                       {return rotate(pos_cm_old, rot_angle_vec, x, rot_angle_vec.length()) + vel_cm*dt;}); //rotating about the old center of mass, then moving forward

#ifndef NDEBUG
       for (const auto &v : model.vertices) {
           if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
               std::cerr << "rigid body got infinite position for one of its points\n";
           }
       }
#endif

       //updating all normals
       std::transform(model.normals.begin(), model.normals.end(), model.normals.begin(),
                      [&](const vec3& x)
                      {return rotate(pos_cm_old, rot_angle_vec, x, rot_angle_vec.length());}); //rotating about the old center of mass, then moving forward

#ifndef NDEBUG
        for (const auto &v : model.normals) {
            if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
                std::cerr << "rigid body got infinite normal for one of its points\n";
            }
        }
#endif

        //updating all velocities
        std::transform(model.velocities.begin(), model.velocities.end(), model.velocities.begin(),
                       [&](const vec3& x)
                       {return vel_cm + cross( (x-pos_cm),w_cm);}); //velocity of cm plus rotational velocity (v=r x w)
                                                //https://en.wikipedia.org/wiki/Angular_velocity

#ifndef NDEBUG
        for (const auto &v : model.velocities) {
            if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
                std::cerr << "rigid body got infinite normal for one of its points\n";
            }
        }
        for (const auto &v : model.mass) {
            if (!std::isfinite(v) ) {
                std::cerr << "rigid body got infinite mass for one of its points\n";
            }
        }
#endif
    }



    void write_pos(const char* file_loc) {
        std::ofstream output(file_loc);
        for (const auto& p : model.vertices) {
            output << p << "\n";
        }
        output.close();
    }
};







#endif //CODE_BODY_HPP
