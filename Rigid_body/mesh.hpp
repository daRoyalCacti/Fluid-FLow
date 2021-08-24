//
// Created by jacob on 24/8/21.
//

#ifndef CODE_MESH_HPP
#define CODE_MESH_HPP

#include <vector>
#include "../MyMath/vec3.hpp"
#include "triangle.hpp"

struct mesh final{
    std::vector<vec3> vertices;
    std::vector<vec3> velocities;
    std::vector<unsigned> indices;
    const std::vector<double> mass;

    mesh() = delete;
    mesh(const std::vector<vec3> &vertices_, const std::vector<double> &mass_) noexcept : vertices(vertices_), indices{}, mass(mass_), velocities{}  {
#ifndef NDEBUG
        if (vertices_.size() != mass_.size()) {
            std::cerr << "vertices and mass must be the same size\n";
        }
#endif
    }
    mesh(const std::vector<vec3> &vertices_, const std::vector<unsigned> &indices_, const std::vector<double> &mass_) noexcept : vertices(vertices_), indices(indices_), mass(mass_), velocities{} {
#ifndef NDEBUG
        if (vertices_.size() != mass_.size()) {
            std::cerr << "vertices and mass must be the same size\n";
        }
        if (indices_.size() < vertices_.size()) {
            std::cerr << "having less indices that vertices doesn't make sense\n";
        }
#endif
    }

    mesh(const std::vector<vec3> &vertices_, const std::vector<unsigned> &indices_, const std::vector<double> &mass_, const std::vector<vec3> &velocities_) noexcept :
            vertices(vertices_), indices(indices_), mass(mass_), velocities(velocities_) {
#ifndef NDEBUG
        if (vertices_.size() != mass_.size()) {
            std::cerr << "vertices and mass must be the same size\n";
        }
        if (indices_.size() < vertices_.size()) {
            std::cerr << "having less indices that vertices doesn't make sense\n";
        }
        if (vertices_.size() != vertices_.size()) {
            std::cerr << "verticies and velocities must be the same size\n";
        }
#endif
    }
};

#endif //CODE_MESH_HPP
