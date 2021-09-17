//
// Created by jacob on 24/8/21.
//

#ifndef CODE_MESH_HPP
#define CODE_MESH_HPP

#include <vector>
#include <numeric>
#include "../MyMath/vec3.hpp"
#include "triangle.hpp"
#include <numeric>


struct bounding_box {
    vec3 min, max;
};


struct mesh final{
    template <class T>
    class no_init_alloc : public std::allocator<T> {
    public:
        using std::allocator<T>::allocator;

        template <class U, class... Args> void construct(U*, Args&&...) {}
    };

    std::vector<vec3> vertices;
    std::vector<vec3, no_init_alloc<vec3>> velocities;
    std::vector<unsigned> indices;
    const std::vector<double> mass;
    std::vector<vec3> normals;
    vec3 v, w;  //velocity of center of mass and angular velocity about center of mass
    bounding_box bounds;

    mesh() = delete;

    mesh(const std::vector<vec3> &vertices_, const std::vector<unsigned> &indices_, const std::vector<double> &mass_, const std::vector<vec3> &normals_, vec3 v_, vec3 w_) noexcept :
    vertices(vertices_), indices(indices_), mass(mass_), velocities{}, normals(normals_), v(v_), w(w_) {
        velocities.resize(vertices_.size());
        //filling the initial velocities based off the cm velocity and angular velocity
        vec3 pos_cm =  std::inner_product(mass.begin(), mass.end(), vertices.begin(), vec3{}) / std::accumulate(mass.begin(), mass.end(), 0.0);
        double minx=std::numeric_limits<double>::max(), miny=std::numeric_limits<double>::max(), minz=std::numeric_limits<double>::max();
        double maxx=std::numeric_limits<double>::min(), maxy=std::numeric_limits<double>::min(), maxz=std::numeric_limits<double>::min();
        for (size_t i = 0; i < velocities.size(); i++) {
            velocities[i] = v_ + cross( (vertices_[i] - pos_cm), w_);

            if (vertices_[i].x() < maxx) {maxx = vertices_[i].x();}
            if (vertices_[i].x() > minx) {minx = vertices_[i].x();}
            if (vertices_[i].y() < maxy) {maxy = vertices_[i].y();}
            if (vertices_[i].y() > miny) {miny = vertices_[i].y();}
            if (vertices_[i].z() < maxz) {maxz = vertices_[i].z();}
            if (vertices_[i].z() > minz) {minz = vertices_[i].z();}
        }
        bounds = bounding_box{ {minx, miny, minz}, {maxx, maxy, maxz} };


#ifndef NDEBUG
        if (vertices_.size() != mass_.size()) {
            std::cerr << "vertices and mass must be the same size\n";
        }
        if (vertices_.size() != normals_.size()) {
            std::cerr << "vertices and normals must be the same size\n";
        }
        if (indices_.size() < vertices_.size()) {
            std::cerr << "having less indices that vertices doesn't make sense\n";
        }
#endif
    }

    [[nodiscard]] const vec3& get_vertice_index( const unsigned i) const noexcept {
        return vertices[indices[i]];
    }

    [[nodiscard]] const vec3& get_normal_index( const unsigned i) const noexcept {
        return normals[indices[i]];
    }

    [[nodiscard]] const vec3& get_velocity_index( const unsigned i) const noexcept {
        return velocities[indices[i]];
    }
};

#endif //CODE_MESH_HPP
