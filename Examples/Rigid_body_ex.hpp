//
// Created by jacob on 8/8/21.
//


#ifndef CODE_RIGID_BODY_EX_HPP
#define CODE_RIGID_BODY_EX_HPP


#include "../Rigid_body/body.hpp"
#include <string>

#include "../MyMath/ops.hpp"

void Rigid_body_ex() {
    std::vector<vec3> pos;
    std::vector<double> mass;
    std::vector<vec3> forces;

    constexpr unsigned N = 10;
    constexpr double dt = 0.01;
    constexpr double end_time = 5;
    constexpr double timesteps = 10*end_time / dt;

    pos.resize(N*N);
    mass.resize(N*N);
    forces.resize(N*N);

    unsigned counter = 0;
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++, counter++) {
            pos[counter] = vec3(i,j,0);
            mass[counter] = 1;
            forces[counter] = vec3{};
        }
    }

    body b(pos, mass);

    unsigned writing_counter = 0;

    forces[0] = vec3(0,0,N*N);
    //forces[N-1] = vec3(static_cast<double>(N*N)/10,0,0);
    forces[N*N-1] = vec3(0,0,-static_cast<double>(N*N));

    for (unsigned i = 0; i < timesteps; i++) {
        /*if (i== 0 ) {
            forces[0] = vec3(0,0,N*N);
            forces[N-1] = vec3(0,0,-static_cast<double>(N*N));
        }
        if (i==1) {
            forces[0] = vec3(0,0,0);
            forces[N-1] = vec3(0,0,0);
        }*/
        b.update_pos(forces, b.model.vertices, dt);

        if (i % static_cast<unsigned >(0.1/dt) == 0) {
            std::string file_name;
            if (writing_counter < 10) {
                file_name = "000" + std::to_string(writing_counter++);
            } else if (writing_counter < 100) {
                file_name = "00" + std::to_string(writing_counter++);
            } else if (writing_counter < 1000) {
                file_name = "0" + std::to_string(writing_counter++);
            } else {
                file_name = std::to_string(writing_counter++);
            }
            b.write_pos(("../DEBUG/rigid_body_data/" + file_name + ".txt").c_str());
        }
    }
}

#endif //CODE_RIGID_BODY_EX_HPP
