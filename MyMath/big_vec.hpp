//
// Created by jacob on 9/8/21.
//

#ifndef CODE_BIG_VEC_HPP
#define CODE_BIG_VEC_HPP

#include <Eigen/Dense>

#include <algorithm>

#include "vec3.hpp"
#include "finite_difference.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include "../Fluid_flow/boundary_conditions.hpp"

template <typename T>
struct big_vec {
    const grid* g;

    big_vec() : g{} {}
    explicit big_vec(const boundary_conditions &b) noexcept  : g{&b.global_grid} {}


    [[nodiscard]] inline bool has_left(const unsigned ind) const noexcept {return g->has_left(ind);}
    [[nodiscard]] inline bool has_right(const unsigned ind) const noexcept {return g->has_right(ind);}
    [[nodiscard]] inline bool has_down(const unsigned ind) const noexcept {return g->has_down(ind);}
    [[nodiscard]] inline bool has_up(const unsigned ind) const noexcept {return g->has_up(ind);}
    [[nodiscard]] inline bool has_front(const unsigned ind) const noexcept {return g->has_front(ind);}
    [[nodiscard]] inline bool has_back(const unsigned ind) const noexcept {return g->has_back(ind);}


    [[nodiscard]] inline bool has_2left(const unsigned ind) const noexcept {return g->has_2left(ind);}
    [[nodiscard]] inline bool has_2right(const unsigned ind) const noexcept {return g->has_2right(ind);}
    [[nodiscard]] inline bool has_2down(const unsigned ind) const noexcept {return g->has_2down(ind);}
    [[nodiscard]] inline bool has_2up(const unsigned ind) const noexcept {return g->has_2up(ind);}
    [[nodiscard]] inline bool has_2front(const unsigned ind) const noexcept {return g->has_2front(ind);}
    [[nodiscard]] inline bool has_2back(const unsigned ind) const noexcept {return g->has_2back(ind);}

    [[nodiscard]] inline bool is_boundary(const unsigned ind) const noexcept { return g->is_boundary(ind); }
    [[nodiscard]] inline bool is_inside_boundary(const unsigned ind) const noexcept { return g->is_inside_boundary(ind); }

    [[nodiscard]] constexpr inline double dx(const unsigned ind) const noexcept {return g->dx;}
    [[nodiscard]] constexpr inline double dy(const unsigned ind) const noexcept {return g->dy;}
    [[nodiscard]] constexpr inline double dz(const unsigned ind) const noexcept {return g->dz;}
    [[nodiscard]] inline vec3 get_pos(const unsigned ind) const noexcept {return (*g)[ind];}
    [[nodiscard]] inline unsigned get_inds(const vec3& p) const noexcept {    std::cerr << "calling get_inds for big vec --- function is very slow\n";
        return g->get_ind(p);
    }


    virtual unsigned long size() const {
#ifndef NDEBUG
        std::cerr << "virtual size function of big_vec called. This should never happen\n";
#endif
        return -1;  //will overflow and give a stupid number
    }
    virtual void clear() {
#ifndef NDEBUG
        std::cerr << "virtual clear function of big_vec called. This should never happen\n";
#endif
    }


    [[nodiscard]] virtual inline T move(const unsigned ind, const int x, const int y, const int z) const noexcept {
#ifndef NDEBUG
        std::cerr << "virtual move function of big_vec called. This should never happen\n";
#endif
        return T{};
    }

    [[nodiscard]] inline T move(const unsigned ind, const vec3&v) const noexcept {
        return move(ind, v.x(), v.y(), v.z());
    }

    [[nodiscard]] virtual T operator()(const unsigned ind) const noexcept {
#ifndef NDEBUG
        std::cerr << "virtual () overload of big_vec called. This should never happen\n";
#endif
        T a{};
        return a;
    }
    [[nodiscard]] virtual T& operator()(const unsigned ind) noexcept {
#ifndef NDEBUG
        std::cerr << "virtual () overload of big_vec called. This should never happen\n";
#endif
        T a{};
        return a;
    }

    [[nodiscard]] unsigned get_move_ind(const unsigned ind, const int x, const int y, const int z) const noexcept {
        auto curr_ind = ind;
#ifndef NDEBUG
        if (curr_ind == -1) {
            std::cerr << "index is -1. This should never happen.\n";
        }
        if (curr_ind >= size()) {
            std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
        }
#endif

        if (x > 0) {
            for (unsigned i = 0; i < x; i++) {
                curr_ind = g->r[curr_ind].right;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "x move failed. current x : " << i+1 << "/" << x << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "x move failed. current x : " << i+1 << "/" << x << "\n";
                }
#endif
            }
        } else {
            for (int i = x; i < 0; i++) {
                curr_ind = g->r[curr_ind].left;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "x move failed. current x : " << i+1 << "/" << x << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "x move failed. current x : " << i+1 << "/" << x << "\n";
                }
#endif
            }
        }
        if (y > 0) {
            for (unsigned j = 0; j < y; j++) {
                curr_ind = g->r[curr_ind].up;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "y move failed. current y : " << j+1 << "/" << y << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "y move failed. current y : " << j+1 << "/" << y << "\n";
                }
#endif
            }
        } else {
            for (int j = y; j < 0; j++) {
                curr_ind = g->r[curr_ind].down;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "y move failed. current y : " << j+1 << "/" << y << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "y move failed. current y : " << j+1 << "/" << y << "\n";
                }
#endif
            }
        }
        if (z > 0) {
            for (unsigned k = 0; k < z; k++) {
                curr_ind = g->r[curr_ind].back;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "z move failed. current z : " << k+1 << "/" << z << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "z move failed. current z : " << k+1 << "/" << z << "\n";
                }
#endif
            }
        } else {
            for (int k = z; k < 0; k++) {
                curr_ind = g->r[curr_ind].front;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "z move failed. current z : " << k+1 << "/" << z << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "z move failed. current z : " << k+1 << "/" << z << "\n";
                }
#endif
            }
        }

        return curr_ind;
    }
};


struct big_vec_d final : public big_vec<double> {
    Eigen::Matrix<double, Eigen::Dynamic, 1> v;


    big_vec_d() : v{}, big_vec() {  }
    explicit big_vec_d(const boundary_conditions &b) noexcept  : big_vec(b) {
        v.resize( static_cast<long>(g->size()) );
        clear(); //ensuring nothing weird happens
    }

    big_vec_d &operator+=(const big_vec_d &other) {
        v += other.v;
        return *this;
    }

    big_vec_d& operator=(const big_vec_d& q) noexcept {
        v = q.v;
        g = q.g;
        return *this;
    }

    [[nodiscard]] unsigned long size() const override {
        return v.size();
    }

    void clear() override {
        for (unsigned i = 0; i < size(); i++) {
            v(i) = 0;
        }
    }



    [[nodiscard]] inline double move(const unsigned ind, const int x, const int y, const int z) const noexcept override {
        const auto new_ind = get_move_ind(ind, x, y, z);

        return v[new_ind];
    }

    [[nodiscard]] double& operator()(const unsigned ind) noexcept override { return v(ind); }
    [[nodiscard]] double operator()(const unsigned ind) const noexcept override { return v(ind); }
};




struct big_vec_v final : public big_vec<vec3> {
    big_vec_d xv, yv, zv;   //needed because these are to be passed into derivative functions


    big_vec_v() : xv{}, yv{}, zv{}, big_vec() {}
    explicit big_vec_v(const boundary_conditions &b) noexcept : big_vec(b), xv(b), yv(b), zv(b) {}

    void clear() override {
        xv.clear();
        yv.clear();
        zv.clear();
    }

    unsigned long size() const override {
        return xv.size();
    }

    big_vec_v& operator=(const big_vec_v& v) noexcept {
        g = v.g;

        xv = v.xv;
        yv = v.yv;
        zv = v.zv;
        return *this;
    }


    [[nodiscard]] inline vec3 move(const unsigned ind, const int x, const int y, const int z) const noexcept override {
        const auto new_ind = get_move_ind(ind, x, y, z);

        return {xv.v[new_ind], yv.v[new_ind], zv.v[new_ind]};
    }

    //not ideal, should make it so () can return a reference
    auto add_elm(const unsigned ind, const vec3 &elm) noexcept {
        return add_elm(ind, elm.x(), elm.y(), elm.z());
    }
    //not ideal, should make it so () can return a reference
    void add_elm(const unsigned ind, const double ex, const double ey, const double ez) noexcept {
        xv(ind) = ex;
        yv(ind) = ey;
        zv(ind) = ez;
    }

    [[nodiscard]] vec3 operator()(const unsigned ind) const noexcept override {
        return {xv.v[ind], yv.v[ind], zv.v[ind]};
    }

};




template<typename T>
void write_vec(const T& v, const char* file_loc) noexcept {
    std::ofstream output(file_loc);
    if (output.is_open()) {
        const auto inds = v.g->get_middle_inds();
        for (const auto ind : inds) {
            const auto g = *v.g;
            output << g[ind].x() << " " << g[ind].y() << " " << v(ind) << "\n";
        }
    } else {
        std::cerr << "failed to open file\n";
    }

    output.close();
}

void write_vec(const auto& v, const auto& inds, const char* file_loc) noexcept {
    std::ofstream output(file_loc);
    if (output.is_open()) {
        for (const auto ind : inds) {
            const grid& g = *v.g;
            const auto pos = g[ind];
            output << pos.x() << " " << pos.y() << " " << v(ind) << "\n";
        }
    } else {
        std::cerr << "failed to open file\n";
    }

    output.close();
}



#endif //CODE_BIG_VEC_HPP
