//
// Created by jacob on 19/8/21.
//

#ifndef CODE_BOUNDARY_HPP
#define CODE_BOUNDARY_HPP

#include <vector>
#include <unordered_map>
#include "../MyMath/vec3.hpp"

/*
template <class T>
class no_init_alloc
        : public std::allocator<T>
{
public:
    using std::allocator<T>::allocator;

    template <class U, class... Args> void construct(U*, Args&&...) {}
};
 */

//only have to set those that are false
struct edge_point final {
    bool has_left = true;
    bool has_right = true;
    bool has_front = true;
    bool has_back = true;
    bool has_up = true;
    bool has_down = true;

    edge_point() = default;
};

/*
struct boundary_points final {
    std::vector<edge_point> v;

    explicit boundary_points(const unsigned total_size) {
        v.resize( total_size);
    }

    [[nodiscard]] auto size() const {
        return v.size();
    }

    void clear() {
        for (auto& v_ : v) {
            v_ = edge_point{};
        }
    }

    [[nodiscard]] edge_point& operator()(const unsigned i, const unsigned j, const unsigned k) noexcept { return v[get_index(i,j,k)]; }
    [[nodiscard]] const edge_point& operator()(const unsigned i, const unsigned j, const unsigned k) const noexcept { return v[get_index(i,j,k)]; }

    [[nodiscard]] constexpr inline bool has_left(const unsigned i, const unsigned j, const unsigned k) const noexcept {return v[get_index(i,j,k)].has_left;}
    [[nodiscard]] constexpr inline bool has_right(const unsigned i, const unsigned j, const unsigned k) const noexcept {return v[get_index(i,j,k)].has_right;}
    [[nodiscard]] constexpr inline bool has_down(const unsigned i, const unsigned j, const unsigned k) const noexcept {return v[get_index(i,j,k)].has_down;}
    [[nodiscard]] constexpr inline bool has_up(const unsigned i, const unsigned j, const unsigned k) const noexcept {return v[get_index(i,j,k)].has_up;}
    [[nodiscard]] constexpr inline bool has_front(const unsigned i, const unsigned j, const unsigned k) const noexcept {return v[get_index(i,j,k)].has_front;}
    [[nodiscard]] constexpr inline bool has_back(const unsigned i, const unsigned j, const unsigned k) const noexcept {return v[get_index(i,j,k)].has_back;}

    [[nodiscard]] constexpr inline bool is_boundary(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        const auto end_x = !has_left(i,j,k) || !has_right(i,j,k);
        const auto end_y = !has_down(i,j,k) || !has_up(i,j,k);
        const auto end_z = !has_front(i,j,k) || !has_back(i,j,k);
        return end_x || end_y || end_z;
    }

    [[nodiscard]] constexpr inline bool is_inside_boundary(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        const auto end_x = !has_left(i,j,k) || !has_right(i,j,k);
        const auto end_y = !has_down(i,j,k) || !has_up(i,j,k);
        const auto end_z = !has_front(i,j,k) || !has_back(i,j,k);
        return end_x && end_y && end_z;
    }


    [[nodiscard]] constexpr inline bool has_2left(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        if (i-1 < 0) return false;
        return v[get_index(i,j,k)].has_left && v[get_index(i-1,j,k)].has_left;
    }
    [[nodiscard]] constexpr inline bool has_2right(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        if (i+1 > N) return false;
        return v[get_index(i,j,k)].has_right && v[get_index(i+1,j,k)].has_right;
    }
    [[nodiscard]] constexpr inline bool has_2down(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        if (j-1 < 0) return false;
        return v[get_index(i,j,k)].has_down && v[get_index(i,j-1,k)].has_down;
    }
    [[nodiscard]] constexpr inline bool has_2up(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        if (j+1 > M) return false;
        return v[get_index(i,j,k)].has_up && v[get_index(i,j+1,k)].has_up;
    }
    [[nodiscard]] constexpr inline bool has_2front(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        if (k-1 < 0) return false;
        return v[get_index(i,j,k)].has_front && v[get_index(i,j,k-1)].has_front;
    }
    [[nodiscard]] constexpr inline bool has_2back(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        if (k+1 > P) return false;
        return v[get_index(i,j,k)].has_back && v[get_index(i,j,k+1)].has_back;
    }

private:
    [[nodiscard]] constexpr inline unsigned get_index(const unsigned i, const unsigned j, const unsigned k) const noexcept {
#ifndef NDEBUG
    if (i < 0) {
        std::cerr << "trying to access i < 0\n";
    }
    if (i > N) {
        std::cerr << "trying to access i > N\n";
    }
    if (j < 0) {
        std::cerr << "trying to access j < 0\n";
    }
    if (j > M) {
        std::cerr << "trying to access j > M\n";
    }
    if (k < 0) {
        std::cerr << "trying to access k < 0\n";
    }
    if (k > P) {
        std::cerr << "trying to access k > P\n";
    }
#endif
        return i + (N+1)*j + (N+1)*(M+1)*k;
    }

};
*/

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
};


struct vel_points final {
    std::unordered_map<unsigned, vec3> m{};   //takes a square index and returns a point on the velocity of the mehs

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

};

#endif //CODE_BOUNDARY_HPP
