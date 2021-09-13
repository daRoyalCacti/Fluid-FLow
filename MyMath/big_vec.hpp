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
    [[nodiscard]] inline vec3 get_pos(const unsigned ind) const noexcept {
        return (*g)[ind];}
        [[nodiscard]] inline unsigned get_inds(const vec3& p) const noexcept {    std::cerr << "calling get_inds for big vec --- function is very slow\n";
        return g->get_ind(p);
    }


    virtual unsigned long size() {
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


    [[nodiscard]] virtual constexpr inline T move(const unsigned ind, const int x, const int y, const int z) const noexcept {
#ifndef NDEBUG
        std::cerr << "virtual move function of big_vec called. This should never happen\n";
#endif
        return T{};
    }

    [[nodiscard]] virtual double& operator()(const unsigned ind) noexcept {
#ifndef NDEBUG
        std::cerr << "virtual () overload of big_vec called. This should never happen\n";
#endif
        double a = 1;
        return a;
    }
    [[nodiscard]] virtual const double& operator()(const unsigned ind) const noexcept {
#ifndef NDEBUG
        std::cerr << "virtual () overload of big_vec called. This should never happen\n";
#endif
        double a = 1;
        return a;
    }

};


struct big_vec_d final : public big_vec<double> {
    Eigen::Matrix<double, Eigen::Dynamic, 1> v;


    big_vec_d() : v{}, big_vec() {  }
    explicit big_vec_d(const boundary_conditions &b) noexcept  : big_vec(b) {
        v.resize( g->size() );
        for (unsigned i = 0; i < v.size(); i++) {
            v(i) = double{};   //ensuring nothing strange
        }
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

    unsigned long size() override {
        return v.size();
    }

    void clear() override {
        for (unsigned i = 0; i < size(); i++) {
            v(i) = 0;
        }
    }



    [[nodiscard]] constexpr inline double move(const unsigned ind, const int x, const int y, const int z) const noexcept override {
        auto curr_ind = ind;

        if (x > 0) {
            for (unsigned i = 0; i < x; i++) {
                curr_ind = g->r[curr_ind].left;
            }
        } else {
            for (unsigned i = -x; i < 0; i++) {
                curr_ind = g->r[curr_ind].right;
            }
        }
        if (y > 0) {
            for (unsigned j = 0; j < y; j++) {
                curr_ind = g->r[curr_ind].up;
            }
        } else {
            for (unsigned j = -y; j < 0; j++) {
                curr_ind = g->r[curr_ind].down;
            }
        }
        if (z > 0) {
            for (unsigned k = 0; k < z; k++) {
                curr_ind = g->r[curr_ind].front;
            }
        } else {
            for (unsigned k = -z; k < 0; k++) {
                curr_ind = g->r[curr_ind].back;
            }
        }

        return v[curr_ind];
    }

    [[nodiscard]] double& operator()(const unsigned ind) noexcept { return v(ind); }
    [[nodiscard]] const double& operator()(const unsigned ind) const noexcept { return v(ind); }
};



/*

struct big_vec_v final : public big_vec {
    big_vec<double> xv, yv, zv;


    big_vec() : dxb{}, dyb{}, dzb{}, xv{}, yv{}, zv{} {}
    big_vec(const double dx_, const double dy_, const double dz_, const boundary_points<N,M,P>* const b_) noexcept : dxb(dx_), dyb(dy_), dzb(dz_),
        xv(dx_, dy_, dz_, b_),
        yv(dx_, dy_, dz_, b_),
        zv(dx_, dy_, dz_, b_) {}

    [[maybe_unused]] big_vec(const big_vec<N, M, P, vec3>& v) noexcept : dxb(v.dx), dyb(v.dy), dzb(v.dz) {
        xv.v = v.xv.v;
        yv.v = v.yv.v;
        zv.v = v.zv.v;
    }

    void clear() {
        xv.clear();
        yv.clear();
        zv.clear();
    }

    big_vec& operator=(const big_vec<N, M, P, vec3>& v) noexcept {
        dxb = v.dxb;
        dyb = v.dyb;
        dzb = v.dzb;

        xv = v.xv;
        yv = v.yv;
        zv = v.zv;
        return *this;
    }

    [[nodiscard]] constexpr inline bool has_left(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_left(i,j,k);}
    [[nodiscard]] constexpr inline bool has_right(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_right(i,j,k);}
    [[nodiscard]] constexpr inline bool has_down(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_down(i,j,k);}
    [[nodiscard]] constexpr inline bool has_up(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_up(i,j,k);}
    [[nodiscard]] constexpr inline bool has_front(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_front(i,j,k);}
    [[nodiscard]] constexpr inline bool has_back(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_back(i,j,k);}


    [[nodiscard]] constexpr inline bool has_2left(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_2left(i,j,k);}
    [[nodiscard]] constexpr inline bool has_2right(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_2right(i,j,k);}
    [[nodiscard]] constexpr inline bool has_2down(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_2down(i,j,k);}
    [[nodiscard]] constexpr inline bool has_2up(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_2up(i,j,k);}
    [[nodiscard]] constexpr inline bool has_2front(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_2front(i,j,k);}
    [[nodiscard]] constexpr inline bool has_2back(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.has_2back(i,j,k);}

    [[nodiscard]] constexpr inline vec3 move(const unsigned i, const unsigned j, const unsigned k, const int x, const int y, const int z) const noexcept {
        return vec3(xv.move(i,j,k,x,y,z), yv.move(i,j,k,x,y,z),zv.move(i,j,k,x,y,z));
    }

    [[nodiscard]] constexpr inline double dx(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.dx(i,j,k);}
    [[nodiscard]] constexpr inline double dy(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.dy(i,j,k);}
    [[nodiscard]] constexpr inline double dz(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.dz(i,j,k);}
    [[nodiscard]] constexpr inline vec3 get_pos(const unsigned i, const unsigned j, const unsigned k) const noexcept {return xv.get_pos(i,j,k);}
    [[nodiscard]] constexpr inline vec3 get_inds(const vec3& v) const noexcept {return xv.get_inds(v);}

    [[nodiscard]] constexpr inline bool is_boundary(const unsigned i, const unsigned j, const unsigned k) const noexcept { return xv.is_boundary(i,j,k); }
    [[nodiscard]] constexpr inline bool is_inside_boundary(const unsigned i, const unsigned j, const unsigned k) const noexcept { return xv.is_inside_boundary(i,j,k); }

    //not ideal, should make it so () can return a reference
    void add_elm(const unsigned i, const unsigned j, const unsigned k, const vec3 elm) noexcept {
        xv(i,j,k) = elm.x();
        yv(i,j,k) = elm.y();
        zv(i,j,k) = elm.z();
    }
    //not ideal, should make it so () can return a reference
    void add_elm(const unsigned i, const unsigned j, const unsigned k, const double ex, const double ey, const double ez) noexcept {
        xv(i,j,k) = ex;
        yv(i,j,k) = ey;
        zv(i,j,k) = ez;
    }

    [[nodiscard]] vec3 operator()(const unsigned i, const unsigned j, const unsigned k) const noexcept {
        const auto in = get_index(i,j,k);
        return vec3(xv.v(in), yv.v(in), zv.v(in));
    } //{ return v( get_index(i,j,k) ); }


private:
    double dxb, dyb, dzb;
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




template<unsigned N, unsigned M, unsigned P, typename T>
void write_vec(const big_vec<N,M,P,T>& v, const char* file_loc) noexcept {
    std::ofstream output(file_loc);
    if (output.is_open()) {
        const auto k = P/2;
        for (unsigned i = 0; i <= N; i++) {
            for (unsigned j = 0; j <= M; j++) {
                const auto pos = v.get_pos(i,j,k);
                output << pos.x() << " " << pos.y() << " " << v(i, j, k) << "\n";
            }
        }
    } else {
        std::cerr << "failed to open file\n";
    }

    output.close();
}
*/


#endif //CODE_BIG_VEC_HPP
