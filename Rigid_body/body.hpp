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

    const double M;   //the mass of the system

    body() = delete;

    explicit body(mesh model_) : model(std::move(model_)), M(std::accumulate(model_.mass.begin(), model_.mass.end(), 0.0)),
            pos_cm{ std::inner_product(model_.mass.begin(), model_.mass.end(), model_.vertices.begin(), vec3{}) / std::accumulate(model_.mass.begin(), model_.mass.end(), 0.0)} {
        for (const auto & v : model.velocities) {
            std::cerr << v << "\n";
        }
    }


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
        model.v += a*dt;
        const auto pos_cm_old = pos_cm;
        pos_cm += model.v*dt;

        //finding the torque
        const vec3 t = std::inner_product(points.begin(), points.end(), forces.begin(), vec3{}, std::plus<>(),
                                          [&](const vec3 &r, const vec3 &F){return cross( (r-pos_cm), F);});


        //moment of inertia
        const c_line r_axis(pos_cm, t);
        const double I = std::inner_product(model.mass.begin(), model.mass.end(), model.vertices.begin(), 0.0, std::plus<>(),
                                            [&](const double &m, const vec3 &r){return m * r_axis.distance(r)*r_axis.distance(r);});

        //angular acceleration
        vec3 alpha;
        //can happen due to round off errors when dealing with very small rotations
        if (!std::isfinite(I)) [[unlikely]] {
            alpha = vec3(0);
#ifndef NDEBUG
            std::cerr << "dealing with the invalid I\n";    //lots of errors get printed when this branch is activated, just saying that everything is all good
#endif
        }  else [[likely]] {    //everything all good
            alpha = t/I;
        }


        //finding angle and angular-velocity pseudo-vector
        model.w += alpha*dt;
        const vec3 rot_angle_vec = model.w*dt;

        std::cerr << "\n" << a << "\t" << model.v << "\t" << model.w << "\n";

#ifndef NDEBUG
        {
            bool err = false;
            if (!std::isfinite(model.v.x()) || !std::isfinite(model.v.y()) || !std::isfinite(model.v.z())) {
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
            /*if (!std::isfinite(I)) {
                std::cerr << "rigid body got infinite moment of inertia\n";
                err = true;
            }*/
            if (!std::isfinite(alpha.x()) || !std::isfinite(alpha.y()) || !std::isfinite(alpha.z())) {
                err = true;
                std::cerr << "rigid body got infinite angular acceleration\n";
            }
            if (!std::isfinite(model.w.x()) || !std::isfinite(model.w.y()) || !std::isfinite(model.w.z())) {
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
                std::cerr << "\tvelocity : " << model.v << "\n";
                std::cerr << "\ttorque : " << t << "\n";
                std::cerr << "\tmoment of inertia : " << I << "\n";
                std::cerr << "\tangular accelerator : " << alpha << "\n";
                std::cerr << "\tangular velocity : " << model.w << "\n";
                std::cerr << "\trotation vector : " << rot_angle_vec << "\n";
            }
        }
#endif



        //updating all positions
        std::transform(model.vertices.begin(), model.vertices.end(), model.vertices.begin(),
                       [&](const vec3& x)
                       {return rotate(pos_cm_old, rot_angle_vec, x, rot_angle_vec.length()) + model.v*dt;}); //rotating about the old center of mass, then moving forward

#ifndef NDEBUG
       for (const auto &v : model.vertices) {
           if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
               std::cerr << "rigid body got infinite position for one of its points\n";
           }
       }
#endif

        model.update_bounding_box();

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
                       {return model.v + cross( (x-pos_cm),model.w);}); //velocity of cm plus rotational velocity (v=r x w)
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
