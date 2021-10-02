//
// Created by jacob on 9/8/21.
//

#ifndef CODE_BIG_VEC_HPP
#define CODE_BIG_VEC_HPP

#include <Eigen/Dense>
//#include <Eigen/SparseLU>
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
                    std::cerr << "\tx move failed. current x : " << i+1 << "/" << x << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tx move failed. current x : " << i+1 << "/" << x << "\n";
                }
#endif
            }
        } else {
            for (int i = x; i < 0; i++) {
                curr_ind = g->r[curr_ind].left;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\tx move failed. current x : " << i+1 << "/" << x << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tx move failed. current x : " << i+1 << "/" << x << "\n";
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
                    std::cerr << "\ty move failed. current y : " << j+1 << "/" << y << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\ty move failed. current y : " << j+1 << "/" << y << "\n";
                }
#endif
            }
        } else {
            for (int j = y; j < 0; j++) {
                curr_ind = g->r[curr_ind].down;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\ty move failed. current y : " << j+1 << "/" << y << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\ty move failed. current y : " << j+1 << "/" << y << "\n";
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
                    std::cerr << "\tz move failed. current z : " << k+1 << "/" << z << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tz move failed. current z : " << k+1 << "/" << z << "\n";
                }
#endif
            }
        } else {
            for (int k = z; k < 0; k++) {
                curr_ind = g->r[curr_ind].front;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\tz move failed. current z : " << k+1 << "/" << z << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tz move failed. current z : " << k+1 << "/" << z << "\n";
                }
#endif
            }
        }

        return curr_ind;
    }

    [[nodiscard]] auto get_move_ind(const unsigned ind, const vec3&  v) const noexcept {
        return get_move_ind(ind, v.x(), v.y(), v.z());
    }


    [[nodiscard]] bool can_move(const unsigned ind, const int x, const int y, const int z) const noexcept {
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
                if (curr_ind == -1) {
                    return false;
                }
            }
        } else {
            for (int i = x; i < 0; i++) {
                curr_ind = g->r[curr_ind].left;
                if (curr_ind == -1) {
                    return false;
                }
            }
        }
        if (y > 0) {
            for (unsigned j = 0; j < y; j++) {
                curr_ind = g->r[curr_ind].up;
                if (curr_ind == -1) {
                    return false;
                }
            }
        } else {
            for (int j = y; j < 0; j++) {
                curr_ind = g->r[curr_ind].down;
                if (curr_ind == -1) {
                    return false;
                }
            }
        }
        if (z > 0) {
            for (unsigned k = 0; k < z; k++) {
                curr_ind = g->r[curr_ind].back;
                if (curr_ind == -1) {
                    return false;
                }
            }
        } else {
            for (int k = z; k < 0; k++) {
                curr_ind = g->r[curr_ind].front;
                if (curr_ind == -1) {
                    return false;
                }
            }
        }

        return true;
    }

    [[nodiscard]] auto can_move(const unsigned ind, const vec3&  v) const noexcept {
        return can_move(ind, v.x(), v.y(), v.z());
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


    //interpolates/extrapolates values based on how much the grid moves
    //https://en.wikipedia.org/wiki/Trilinear_interpolation#Alternative_algorithm
    // would be quite simple to update to quadratic interpolation
    void move(const double x_off, const double y_off, const double z_off) override {
        //variable to hold the new interpolated values
        Eigen::Matrix<double, Eigen::Dynamic, 1> buffer;
        buffer.resize(static_cast<long>(g->size()));





        //finding points to use for interpolation
        // - assumes the offsets are smaller than the step size
#ifndef NDEBUG
        if (x_off > g->dx) {
            std::cerr << "movement in the x-direction is larger than step size. Interpolation of values might be inaccurate\n";
        }
        if (y_off > g->dy) {
            std::cerr << "movement in the y-direction is larger than step size. Interpolation of values might be inaccurate\n";
        }
        if (z_off > g->dz) {
            std::cerr << "movement in the z-direction is larger than step size. Interpolation of values might be inaccurate\n";
        }
#endif
        constexpr unsigned no_points = 20;
        for (unsigned index = 0; index < buffer.size(); index++) {
            //std::cerr << index << "/" << buffer.size() -1 << "\n";

            unsigned interp_indices[no_points];
            size_t counter = 0;
            /*size_t counter = 1;

            interp_indices[0] = index;  //use itself in innterpolation

            //setting the indices along the major axes
            for (const auto& horiz : {-1,1}) {
                if (can_move(index, horiz, 0, 0)) {
                    interp_indices[counter++] = get_move_ind(index, horiz, 0, 0);
                }
            }
            for (const auto& vert : {-1,1}) {
                if (can_move(index, 0, vert, 0)) {
                    interp_indices[counter++] = get_move_ind(index, 0, vert, 0);
                }
            }
            for (const auto& in : {-1,1}) {
                if (can_move(index, 0, 0, in)) {
                    interp_indices[counter++] = get_move_ind(index, 0, 0, in);
                }
            }


            //getting end values
            for (const auto& horiz : {-1,0,1}) {
                for (const auto& vert : {-1,1}) {
                    if (can_move(index, horiz, vert, 0)) {
                        interp_indices[counter++] = get_move_ind(index, horiz, vert, 0);
                        if (counter == 8) {
                            goto got_indices;
                        }
                    }
                }
            }

            for (const auto& horiz : {-1,0,1}) {
                for (const auto& in : {-1,1}) {
                    if (can_move(index, horiz, 0, in)) {
                        interp_indices[counter++] = get_move_ind(index, horiz, 0, in);
                        if (counter == 8) {
                            goto got_indices;
                        }
                    }
                }
            }

            //getting corner values
            for (const auto& horiz : {-1,0,1}) {
                for (const auto& vert : {-1,1}) {
                    for (const auto& in : {-1,1}) {
                        if (can_move(index, horiz, vert, in)) {
                            interp_indices[counter++] = get_move_ind(index, horiz, vert, in);
                            if (counter == 8) {
                                goto got_indices;
                            }
                        }
                    }
                }
            }



            got_indices:
#ifndef NDEBUG
            if (counter != 8) {
                std::cerr << "move points are needed for interpolation\n";
            }
#endif
             */
            int i = 1;
            while (true) {
                //std::cerr << "\t" << counter << "\n";
                //getting corner values first
                for (const auto& horiz : {-i,i}) {
                    for (const auto& vert : {-i,i}) {
                        for (const auto& in : {-i,i}) {
                            if (can_move(index, horiz, vert, in)) {
                                interp_indices[counter++] = get_move_ind(index, horiz, vert, in);
                                if (counter == no_points) {
                                    goto got_indices;
                                }
                            }
                        }
                    }
                }

                //then getting edge values
                for (const auto& horiz : {-i,i}) {
                    for (const auto& vert : {-i,i}) {
                        if (can_move(index, horiz, vert, 0)) {
                            interp_indices[counter++] = get_move_ind(index, horiz, vert, 0);
                            if (counter == no_points) {
                                goto got_indices;
                            }
                        }
                    }

                    for (const auto& in : {-i,i}) {
                        if (can_move(index, horiz, 0, in)) {
                            interp_indices[counter++] = get_move_ind(index, horiz, 0, in);
                            if (counter == no_points) {
                                goto got_indices;
                            }
                        }
                    }
                }


                //then points along major axes
                for (const auto& horiz : {-i,i}) {
                    if (can_move(index, horiz, 0, 0)) {
                        interp_indices[counter++] = get_move_ind(index, horiz, 0, 0);
                        if (counter == no_points) {
                            goto got_indices;
                        }
                    }
                }
                for (const auto& vert : {-i,i}) {
                    if (can_move(index, 0, vert, 0)) {
                        interp_indices[counter++] = get_move_ind(index, 0, vert, 0);
                        if (counter == no_points) {
                            goto got_indices;
                        }
                    }
                }
                for (const auto& in : {-i,i}) {
                    if (can_move(index, 0, 0, in)) {
                        interp_indices[counter++] = get_move_ind(index, 0, 0, in);
                        if (counter == no_points) {
                            goto got_indices;
                        }
                    }
                }



                i++;
            }
            got_indices:





            //finding the constants in y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
            //setting the matrix
            //Eigen::Matrix<double, 8, 8> mat;
            Eigen::Matrix<double, no_points, 1> vec;
            Eigen::SparseMatrix<double> mat(no_points,no_points);
            //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

            for (unsigned i = 0; i < no_points; i++) {
                //std::cerr << i << "\n";
                vec[i] = v[interp_indices[i]];

                const auto x = g->x[interp_indices[i]];
                const auto y = g->y[interp_indices[i]];
                const auto z = g->z[interp_indices[i]];

                mat.insert(i, 0) = 1;
                mat.insert(i, 1) = x;
                mat.insert(i, 2) = y;
                mat.insert(i, 3) = z;
                mat.insert(i, 4) = x*y;
                mat.insert(i, 5) = x*z;
                mat.insert(i, 6) = y*z;
                mat.insert(i, 7) = x*y*z;

                mat.insert(i, 8) = x*x;
                mat.insert(i, 9) = y*y;
                mat.insert(i, 10) = z*z;
                mat.insert(i, 11) = x*x*y;
                mat.insert(i, 12) = x*x*z;
                mat.insert(i, 13) = y*y*x;
                mat.insert(i, 14) = y*y*z;
                mat.insert(i, 15) = x*z*z;
                mat.insert(i, 16) = y*z*z;
                mat.insert(i, 17) = x*x*y*z;
                mat.insert(i, 18) = x*y*y*z;
                mat.insert(i, 19) = x*y*z*x;
            }

            /*if (index == 565) {//43) {
                for (unsigned i = 0; i < 8; i++) {
                    std::cerr << v[interp_indices[i]] << "\t";
                }
                std::cerr << "\n";

                for (unsigned i = 0; i < 8; i++) {
                    std::cerr << g->x[interp_indices[i]] << "\t";
                }
                std::cerr << "\t (" << g->x[index] + x_off << ")\n";

                for (unsigned i = 0; i < 8; i++) {
                    std::cerr << g->y[interp_indices[i]] << "\t";
                }
                std::cerr << "\t (" << g->y[index] + y_off << ")\n";

                for (unsigned i = 0; i < 8; i++) {
                    std::cerr << g->z[interp_indices[i]] << "\t";
                }
                std::cerr << "\t (" << g->z[index] + z_off << ")\n";
                std::cerr << "\n";
            }*/
            //mat and vec now set, just need to solve for the coefficients
            //solver.analyzePattern(mat);
            //solver.factorize(mat);
            //std::cerr << 1 << "\n";
            solver.compute(mat);
            //std::cerr << 2 << "\n";
            const Eigen::Matrix<double, no_points, 1> a = solver.solve(vec);

            //const Eigen::Matrix<double, 8, 1> a = mat.inverse()*vec;

            const auto x = g->x[index] + x_off;
            const auto y = g->y[index] + y_off;
            const auto z = g->z[index] + z_off;
            //std::cerr << 3 << "\n";
            buffer[index] = a(0) + a(1)*x + a(2)*y + a(3)*z + a(4)*x*y + a(5)*x*z + a(6)*y*z + a(7)*x*y*z +
                    a(8)*x*x + a(9)*y*y + a(10)*z*z +  a(11)*x*x*y + a(12)*x*x*z + a(13)*y*y*x +
                    a(14)*y*y*z + a(15)*x*z*z + a(16)*y*z*z + a(17)*x*x*y*z + a(18)*x*y*y*z + a(19)*x*y*z*x ;

            //std::cerr << 4 << "\n";
        }

        v = std::move(buffer);

    }
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


    void move(const double x_off, const double y_off, const double z_off) override {
        xv.move(x_off, y_off, z_off);
        yv.move(x_off, y_off, z_off);
        zv.move(x_off, y_off, z_off);
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

/*
//https://en.wikipedia.org/wiki/Trilinear_interpolation#Alternative_algorithm
// uses linear interpolation to update points that were previously inside a boundary
// would be quite simple to update to quadratic interpolation
// uses interpolation instead of an explicit extrapolation procedure is because I can't find any extrapolation procedures
template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::extrapolate(big_vec<N,M,P, double> &p) {
    //finds the constants in
    //y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz

    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <=P; k++) {
                //const bool is_bound = bound.is_boundary(i,j,k);
                //const bool was_bound = norms.contains(i,j,k) && norms.normal(i,j,k) != vec3(0); //bound_prev.is_boundary(i,j,k);

                const bool is_bound = norms.contains(i,j,k);
                const bool was_bound = norms_prev.contains(i,j,k) && norms_prev.normal(i,j,k) == vec3(0);


                if (!is_bound && was_bound ) {   //was inside mesh but is current outside of the boundary
                    //the points where the extrapolation is to be used
                    const unsigned xi = i;
                    const unsigned yi = j;
                    const unsigned zi = k;





                    Eigen::Matrix<double, 8, 8> mat;
                    Eigen::Matrix<double, 8, 1> vec;

                    //using values on the faces of ever-increasing cubes
                    for (int s = 2; s < 80; s+=2) {    //upper bound on s arbitrarily large
                        std::vector<unsigned> vals_b, vals_s;
                        vals_b.resize(s+1);
                        vals_s.resize(s-1);
                        std::iota(vals_b.begin(), vals_b.end(), -s/2);  //include the entire surface
                        std::iota(vals_s.begin(), vals_s.end(), -s/2+1);//surface less 1 edge

                        unsigned counter = 0;

                        //moving along the surfaces of the cubes
                        {
                            int x = -s/2;
                            for (const auto &y : vals_b) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int x = s/2;
                            for (const auto &y : vals_b) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int y = -s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int y = s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int z = -s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &y : vals_s) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int z = s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &y : vals_s) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }

                        if (counter >= 8) {
                            break;
                        }


                    }   //end s for loop

                    //mat and vec now set, just need to solve for the coefficients
                    const Eigen::Matrix<double, 8, 1> a = mat.inverse()*vec;

                    //setting the value of p
                    //y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
                    const auto pos = vel_bc.get_pos(xi,yi,zi);
                    const auto x = pos.x();
                    const auto y = pos.y();
                    const auto z = pos.z();
                    p(xi, yi, zi) = a(0) + a(1)*x + a(2)*y + a(3)*z + a(4)*x*y + a(5)*x*z + a(6)*y*z + a(7)*x*y*z;

                }

            }
        }
    }

}

template<typename T>
void big_vec<T>::set_matrix_row(const unsigned x, const unsigned y, const unsigned z, unsigned &counter, Eigen::Matrix<double, 8, 8> &mat,
                                                Eigen::Matrix<double, 8, 1> &vec, const big_vec<T> &p, const double xi, const double yi, const double zi) {
    const auto xp = x+xi;
    const auto yp = y+yi;
    const auto zp = z+zi;

    std::cerr << "set matrix row will soon be wrong\n";


    //https://en.wikipedia.org/wiki/Trilinear_interpolation#Alternative_algorithm
    bool has_normal = norms.contains(xp, yp, zp);
    if (!has_normal || (norms.normal(xp, yp, zp) != vec3(0))  ) {    //if point not inside a boundary
        const auto pos = vel_bc.get_pos(xp, yp, zp);
        const auto x0 = pos.x();
        const auto y0 = pos.y();
        const auto z0 = pos.z();

        mat(counter, 0) = 1;
        mat(counter, 1) = x0;
        mat(counter, 2) = y0;
        mat(counter, 3) = z0;
        mat(counter, 4) = x0*y0;
        mat(counter, 5) = x0*z0;
        mat(counter, 6) = y0*z0;
        mat(counter, 7) = x0*y0*z0;

        vec(counter) = p(xp, yp, zp);
    }

    counter++;
}
*/

#endif //CODE_BIG_VEC_HPP
