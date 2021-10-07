//
// Created by jacob on 24/8/21.
//

#ifndef CODE_MESH_HPP
#define CODE_MESH_HPP

#include <utility>
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

    mesh(const std::vector<vec3> &vertices_, std::vector<unsigned> indices_, std::vector<double> mass_, std::vector<vec3> normals_, vec3 v_, vec3 w_) noexcept :
    vertices(vertices_), indices(std::move(indices_)), mass(std::move(mass_)), velocities{}, normals(std::move(normals_)), v(v_), w(w_) {
        velocities.resize(vertices_.size());
        //filling the initial velocities based off the cm velocity and angular velocity
        vec3 pos_cm =  std::inner_product(mass.begin(), mass.end(), vertices.begin(), vec3{}) / std::accumulate(mass.begin(), mass.end(), 0.0);
        for (size_t i = 0; i < velocities.size(); i++) {
            velocities[i] = v_ + cross( (vertices_[i] - pos_cm), w_);
        }
        update_bounding_box();


#ifndef NDEBUG
        if (vertices.size() != mass.size()) {
            std::cerr << "vertices and mass must be the same size\n";
        }
        if (vertices.size() != normals.size()) {
            std::cerr << "vertices and normals must be the same size\n";
        }
        if (indices.size() < vertices.size()) {
            std::cerr << "having less indices that vertices doesn't make sense\n";
        }
#endif
    }

    mesh(const mesh &m) : mass(m.mass), vertices(m.vertices), indices(m.indices), normals(m.normals), v(m.v), w(m.w), bounds(m.bounds)  {
        velocities.resize(m.velocities.size() );
        for (size_t i = 0; i < m.velocities.size(); i++) {
            velocities[i] = m.velocities[i];
        }

    }

    void update_bounding_box() {
        double minx=std::numeric_limits<double>::max(), miny=std::numeric_limits<double>::max(), minz=std::numeric_limits<double>::max();
        double maxx=std::numeric_limits<double>::min(), maxy=std::numeric_limits<double>::min(), maxz=std::numeric_limits<double>::min();
        for (auto &vert : vertices) {
            if (vert.x() > maxx) {maxx = vert.x();}
            if (vert.x() < minx) {minx = vert.x();}
            if (vert.y() > maxy) {maxy = vert.y();}
            if (vert.y() < miny) {miny = vert.y();}
            if (vert.z() > maxz) {maxz = vert.z();}
            if (vert.z() < minz) {minz = vert.z();}
        }
        bounds = bounding_box{ {minx, miny, minz}, {maxx, maxy, maxz} };
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
