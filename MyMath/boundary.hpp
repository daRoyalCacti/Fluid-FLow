//
// Created by jacob on 19/8/21.
//

#ifndef CODE_BOUNDARY_HPP
#define CODE_BOUNDARY_HPP

#include <vector>
#include <unordered_map>
#include "../MyMath/vec3.hpp"
#include "../MyMath/ops.hpp"
#include "../Rigid_body/triangle.hpp"
#include "../Rigid_body/triangle_mesh.hpp"

struct tri_inds final {
    std::unordered_map<unsigned, const triangle*> m{};   //takes a square index and returns a point on the boundary of a mesh

    tri_inds() = default;
    explicit tri_inds(size_t num_points) noexcept {
        m.reserve(num_points);
    }

    auto size() const {
        return m.size();
    }

    void add_point(const unsigned index, const triangle* t) noexcept {
        m.insert({index, t});
    }

    void clear() {
        m.clear();
    }

    [[nodiscard]] bool contains(const unsigned index) const noexcept {
        return m.contains(index);
    }

    [[nodiscard]] auto get_index(const unsigned index) const noexcept {
        return m.at(index);
    }



};


struct mesh_points final {
    std::unordered_map<unsigned, vec3> m{};   //takes a square index and returns a point on the boundary of a mesh

    mesh_points() = default;
    explicit mesh_points(size_t num_points) noexcept {
        m.reserve(num_points);
    }

    auto size() const {
        return m.size();
    }

    void add_point(const unsigned index, const vec3& normal) noexcept {
        m.insert({index, normal});
    }

    void clear() {
        m.clear();
    }

    [[nodiscard]] bool contains(const unsigned index) const noexcept {
        return m.contains(index);
    }

    [[nodiscard]] vec3 get_point(const unsigned index) const noexcept {
        return m.at(index);
    }

    void update(const vec3 &w, const vec3 &v, const vec3& c_o_m, const double dt) {
        const auto rot_angle_vec = w*dt;
        /*std::transform(model.vertices.begin(), model.vertices.end(), model.vertices.begin(),
                       [&](const vec3& x)
                       {return rotate(pos_cm_old, rot_angle_vec, x, rot_angle_vec.length()) + model.v*dt;}); //rotating about the old center of mass, then moving forward.*/
        for (auto & V : m) {
            V.second = rotate(c_o_m, rot_angle_vec, V.second, rot_angle_vec.length()) + v*dt;
        }
    }


};



struct boundary_normals final {
    std::unordered_map<unsigned, vec3> m{};   //takes a square index and returns the normal vector

    boundary_normals() = default;
    explicit boundary_normals(size_t num_points) noexcept {
        m.reserve(num_points);
    }

    auto size() const {
        return m.size();
    }

    void add_point(const unsigned index, const vec3& normal) noexcept {
        m.insert({index, normal});
    }

    void clear() {
        m.clear();
    }

    [[nodiscard]] bool contains(const unsigned index) const noexcept {
        return m.contains(index);
    }

    [[nodiscard]] vec3 normal(const unsigned index) const noexcept {
        return m.at(index);
    }


};


struct wall_normals final {
    std::unordered_map<unsigned, vec3> m{};   //takes a square index and returns the normal vector

    wall_normals() = default;
    explicit wall_normals(size_t num_points) noexcept {
        m.reserve(num_points);
    }

    auto size() const {
        return m.size();
    }

    void add_point(const unsigned index, const vec3& normal) noexcept {
        m.insert({index, normal});
    }

    void clear() {
        m.clear();
    }

    [[nodiscard]] bool contains(const unsigned index) const noexcept {
        return m.contains(index);
    }

    [[nodiscard]] vec3 normal(const unsigned index) const noexcept {
        return m.at(index);
    }

};


struct boundary_normals_grouped final {
    boundary_normals &b1;
    wall_normals &b2;

    boundary_normals_grouped() = delete;
    boundary_normals_grouped(boundary_normals& b1_,  wall_normals& b2_) : b1(b1_), b2(b2_) {}

    [[nodiscard]] auto size() const {
        return b1.size() + b2.size();
    }


    [[nodiscard]] bool contains(const unsigned index) const noexcept {
        return b1.contains(index) || b2.contains(index);
    }

    [[nodiscard]] vec3 normal(const unsigned index) const noexcept {
        if (b1.contains(index)) {
            return b1.normal(index);
        } else {
            return b2.normal(index);
        }
    }
    
};

struct vel_points final {
    std::unordered_map<unsigned, vec3> m{};   //takes a square index and returns a point on the velocity of the mesh

    vel_points() = default;
    explicit vel_points(size_t num_points) noexcept {
        m.reserve(num_points);
    }

    auto size() const {
        return m.size();
    }

    void add_point(const unsigned index, const vec3& normal) noexcept {
        m.insert({index, normal});
    }

    void clear() {
        m.clear();
    }

    [[nodiscard]] bool contains(const unsigned index) const noexcept {
        return m.contains(index);
    }

    [[nodiscard]] vec3 get_vel(const unsigned index) const noexcept {
        return m.at(index);
    }

    void update(const tri_inds &t_inds, const mesh_points &m_points) {
        //t.get_velocity(r.at(time))
        for (auto & v : m) {
            v.second = t_inds.get_index(v.first)->get_velocity(m_points.get_point(v.first));
        }
    }


};







#endif //CODE_BOUNDARY_HPP
