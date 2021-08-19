//
// Created by jacob on 31/7/21.
//

#ifndef CODE_BODY_HPP
#define CODE_BODY_HPP

#include <vector>
#include <memory>
#include "../MyMath/vec3.hpp"

struct body {
    std::shared_ptr<std::vector<vec3>> pos;    //the position of every particle
    const std::vector<double> mass;    //the mass of every particle

    vec3 pos_cm;  //the position of the center of mass
    vec3 vel_cm{};  //the velocity of the center of mass

    const double M;   //the mass of the system

    vec3 w_cm{};    //the angular velocity of the system about the center of mass

    body() = delete;
    body(std::shared_ptr<std::vector<vec3>> &pos_, const std::vector<double>& mass);

    //finds the positions of all the particles undergoing forces <forces> after time dt
    void update_pos(const std::vector<vec3> &forces, double dt);

    void debug_write_pos(const char* file_loc);
};


#endif //CODE_BODY_HPP
