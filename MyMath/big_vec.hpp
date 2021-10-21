//
// Created by jacob on 9/8/21.
//

#ifndef CODE_BIG_VEC_HPP
#define CODE_BIG_VEC_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>

#include "vec3.hpp"
#include "finite_difference.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include "../Fluid_flow/boundary_conditions.hpp"

template <typename T>
struct big_vec {
    grid* g;

    big_vec() : g{} {}
    explicit big_vec(boundary_conditions &b) noexcept  : g{&b.global_grid} {}

    virtual void move(const double x_off, const double y_off, const double z_off) {
        #ifndef NDEBUG
            std::cerr << "virtual move function of big_vec called. This should never happen\n";
        #endif
    }

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


    [[nodiscard]] virtual unsigned long size() const {
#ifndef NDEBUG
        std::cerr << "virtual size function of big_vec called. This should never happen\n";
#endif
        return -1;  //will overflow and give a stupid number
    }
    [[maybe_unused]] virtual void clear() {
#ifndef NDEBUG
        std::cerr << "virtual clear function of big_vec called. This should never happen\n";
#endif
    }

    [[nodiscard]] auto get_move_ind(const unsigned ind, const int x, const int y, const int z) const noexcept {
        return g->get_move_ind(ind, x, y, z);
    }

    [[nodiscard]] auto get_move_ind(const unsigned ind, const vec3&  v) const noexcept {
        return g->get_move_ind(ind, v);
    }


    [[nodiscard]] auto can_move(const unsigned ind, const int x, const int y, const int z) const noexcept {
        return g->can_move(ind, x, y, z);
    }

    [[nodiscard]] auto can_move(const unsigned ind, const vec3&  v) const noexcept {
        return g->can_move(ind, v);
    }


};


struct big_vec_d final : public big_vec<double> {
    Eigen::Matrix<double, Eigen::Dynamic, 1> v;


    big_vec_d() : v{}, big_vec() {  }
    explicit big_vec_d(boundary_conditions &b) noexcept  : big_vec(b) {
        v.resize( static_cast<long>(g->size()) );
        clear(); //ensuring nothing weird happens
    }

    big_vec_d &operator+=(const big_vec_d &other) {
        v += other.v;
        return *this;
    }

    big_vec_d &operator-=(const big_vec_d &other) {
        v -= other.v;
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



    [[nodiscard]] inline double move(const unsigned ind, const int x, const int y, const int z) const noexcept {
        const auto new_ind = get_move_ind(ind, x, y, z);

        return v[new_ind];
    }

    [[nodiscard]] inline double move(const unsigned ind, const vec3& v) const noexcept {
        const auto new_ind = get_move_ind(ind, v);

        return v[new_ind];
    }

    [[nodiscard]] double& operator()(const unsigned ind) noexcept { return v(ind); }
    [[nodiscard]] double operator()(const unsigned ind) const noexcept { return v(ind); }

};




struct big_vec_v final : public big_vec<vec3> {
    big_vec_d xv, yv, zv;   //needed because these are to be passed into derivative functions


    big_vec_v() : xv{}, yv{}, zv{}, big_vec() {}
    explicit big_vec_v(boundary_conditions &b) noexcept : big_vec(b), xv(b), yv(b), zv(b) {}

    void clear() override {
        xv.clear();
        yv.clear();
        zv.clear();
    }

    [[nodiscard]] unsigned long size() const override {
        return xv.size();
    }

    big_vec_v& operator=(const big_vec_v& v) noexcept {
        g = v.g;

        xv = v.xv;
        yv = v.yv;
        zv = v.zv;
        return *this;
    }


    [[nodiscard]] inline vec3 move(const unsigned ind, const int x, const int y, const int z) const noexcept {
        const auto new_ind = get_move_ind(ind, x, y, z);

        return {xv.v[new_ind], yv.v[new_ind], zv.v[new_ind]};
    }

    [[nodiscard]] inline vec3 move(const unsigned ind, const vec3& v) const noexcept {
        const auto new_ind = get_move_ind(ind, v);

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

    [[nodiscard]] vec3 operator()(const unsigned ind) const noexcept {
        return {xv.v[ind], yv.v[ind], zv.v[ind]};
    }

    big_vec_v &operator+=(const big_vec_v &other) {
        xv += other.xv;
        yv += other.yv;
        zv += other.zv;
        return *this;
    }


    big_vec_v &operator-=(const big_vec_v &other) {
        xv -= other.xv;
        yv -= other.yv;
        zv -= other.zv;
        return *this;
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

void write_vec(const big_vec_d& v, const auto& inds, const char* file_loc) noexcept {
    std::ofstream output(file_loc);
    if (output.is_open()) {
        for (const auto ind : inds) {
            const grid& g = *v.g;
            const auto pos = g.get_plot_pos(ind);
            output << pos.x() << " " << pos.y() << " " << v(ind) << "\n";
        }
    } else {
        std::cerr << "failed to open file\n";
    }

    output.close();
}

void write_vec(const big_vec_v& v, const auto& inds, const char* file_loc) noexcept {
    std::ofstream output(file_loc);
    if (output.is_open()) {
        for (const auto ind : inds) {
            const grid& g = *v.g;
            output << g.x[ind] << " " << g.y[ind] << " " << v(ind) << "\n";
        }
    } else {
        std::cerr << "failed to open file\n";
    }

    output.close();
}


#endif //CODE_BIG_VEC_HPP
