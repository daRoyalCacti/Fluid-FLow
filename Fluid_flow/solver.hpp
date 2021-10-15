//
// Created by jacob on 15/8/21.
//

#ifndef CODE_SOLVER_HPP
#define CODE_SOLVER_HPP

#define CL_TARGET_OPENCL_VERSION 300    //the lastest version supported atm
#define VIENNACL_WITH_EIGEN 1
#include <viennacl/vector.hpp>
//#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/gmres.hpp>

#include "../MyMath/big_matrix.hpp"
#include "../MyMath/big_vec.hpp"
#include "update_vecs.hpp"

/*
template <unsigned N, unsigned M, unsigned P>
struct solver {
    viennacl::compressed_matrix<double>  vcl_sparse_matrix;
    //viennacl::vector<double> vcl_rhs;
    viennacl::linalg::gmres_tag my_gmres_tag;

    solver() = delete;
    explicit solver(const big_matrix<N,M,P> &A) {
        vcl_sparse_matrix.resize( (N+1)*(M+1)*(P+1), (N+1)*(M+1)*(P+1) );
        //vcl_rhs.resize( (N+1)*(M+1)*(P+1) );
        viennacl::copy(A.m, vcl_sparse_matrix);

        //setting up the solver
        //could probably mess with the parameters
        my_gmres_tag = viennacl::linalg::gmres_tag(1e-5, 100, 20);// up to 100 iterations, restart after 20 iterations
    }

    void solve(const big_vec<N,M,P,double> &b, big_vec<N,M,P,double> &x) {
        viennacl::vector<double> vcl_rhs( (N+1)*(M+1)*(P+1) );
        viennacl::copy(b.v, vcl_rhs);   //upload vector to device

        //solving
        viennacl::vector<double> res = viennacl::linalg::solve(vcl_sparse_matrix, vcl_rhs, my_gmres_tag);

        //copy result into x
        viennacl::copy(res, x.v);
    }
};
 */

/*
struct solver_vals {
    unsigned max_iterations;
    unsigned dim;
};

#ifdef VIENNACL_WITH_OPENCL
constexpr solver_vals normal_solve{1000, 20};
constexpr solver_vals detailed_solve{5000, 20};
constexpr solver_vals very_detailed_solve{10000, 25};
#endif

#ifdef VIENNACL_WITH_OPENMP
    constexpr solver_vals normal_solve{1000, 1000};
    constexpr solver_vals detailed_solve{5000, 1000};
    constexpr solver_vals very_detailed_solve{10000, 1000};
#endif*/

constexpr unsigned no_solver_choices = 5;
constexpr unsigned max_its_v[no_solver_choices] = {1000, 5000, 10000};
constexpr unsigned dims_v[no_solver_choices]    = {20, 20, 25};

constexpr unsigned max_its_p[no_solver_choices] = {1000, 5000, 10000};
constexpr unsigned dims_p[no_solver_choices]    = {20, 20, 25};


//solves Ax=b for x
//template <unsigned it, unsigned dim>
void solve(const big_matrix &A, const big_vec_d &b, big_vec_d &x, const unsigned it, const unsigned dim) noexcept {
    viennacl::compressed_matrix<double>  vcl_sparse_matrix( b.size(), b.size() );
    //viennacl::coordinate_matrix<double>  vcl_sparse_matrix( (N+1)*(M+1)*(P+1), (N+1)*(M+1)*(P+1) );
    viennacl::vector<double> vcl_rhs( b.size() );

    //copying data into viennacl
    viennacl::copy(A.m, vcl_sparse_matrix);
    viennacl::copy(b.v, vcl_rhs);

    //setting up the solver
    //could probably mess with the parameters
    //viennacl::linalg::gmres_tag my_gmres_tag(1e-5, 100, 20); // up to 100 iterations, restart after 20 iterations
#ifdef VIENNACL_WITH_OPENCL
    viennacl::linalg::gmres_tag my_gmres_tag(1e-100, it, dim);
#endif
#ifdef VIENNACL_WITH_OPENMP
    viennacl::linalg::gmres_tag my_gmres_tag(1e-5, it, dim);
#endif

    /*
#ifdef ACCURATE_SOLVER
    viennacl::linalg::gmres_tag my_gmres_tag(1e-100, 10000, 25); // up to 100 iterations, restart after 20 iterations
#else
    viennacl::linalg::gmres_tag my_gmres_tag(1e-100, 1000, 20);
#endif
     */

    //viennacl::linalg::gmres_tag my_gmres_tag(1e-50, 1000, 20); // up to 100 iterations, restart after 20 iterations

    //solving
    viennacl::vector<double> res = viennacl::linalg::solve(vcl_sparse_matrix, vcl_rhs, my_gmres_tag);

    //copy result into x
    viennacl::copy(res, x.v);
}


void solve_pressure(const boundary_conditions& BC, const big_matrix &Q, const big_vec_d &s, big_vec_d &p_c, big_vec_d &p, const double accuracy_percent) {
    for (unsigned i; i < no_solver_choices; i++) {
        solve(Q, s, p_c, max_its_p[i], dims_p[i]);
        bool accurate = update_pressure_BC(BC, p_c, accuracy_percent);
        if (accurate) {
            p += p_c;
            const bool accurate2 = update_pressure_BC(BC, p, accuracy_percent);
            if (accurate2) {
                std::cerr << "pressure solved using solver" << i << "\n";
                return;
            } else {
                p -= p_c;
            }
        }
    }

    p += p_c;
    std::cerr << "pressure is still not accurate even with most accurate solver\n";

    /*solve<normal_solve.max_iterations, normal_solve.dim>(Q, s, p_c);
    //enforce_PBC(p_c, BC.norms);

    bool accurate = update_pressure_BC(BC, p_c, accuracy_percent);     //TESTING
     if (accurate) {
         p += p_c;
         const bool accurate2 = update_pressure_BC(BC, p, accuracy_percent);
         if (accurate2) {
             return;
         } else {
             p -= p_c;
         }
     }
     solve<detailed_solve.max_iterations, detailed_solve.dim>(Q, s, p_c);
     accurate = update_pressure_BC(BC, p_c, accuracy_percent);
     if (accurate) {
         p += p_c;
         const bool accurate2 = update_pressure_BC(BC, p, accuracy_percent);
         if (accurate2) {
             return;
         } else {
             p -= p_c;
         }
     }
     solve<very_detailed_solve.max_iterations, very_detailed_solve.dim>(Q, s, p_c);
     accurate = update_pressure_BC(BC, p_c, accuracy_percent);
     if (accurate) {
         p += p_c;
         const bool accurate2 = update_pressure_BC(BC, p, accuracy_percent);
         if (accurate2) {
             return;
         } else {
             std::cerr << "pressure is still not accurate even with most accurate solver\n";
         }
     }*/


}

void solve_velocity(const boundary_conditions& BC, const big_matrix &A, const big_vec_v &b, big_vec_v &vc, big_vec_v &v_n, const double accuracy_percent) {
    for (unsigned i; i < no_solver_choices; i++) {
        //solve(Q, s, p_c, max_its[i], dims[i]);
        solve(A, b.xv, vc.xv, max_its_v[i], dims_v[i]);
        solve(A, b.yv, vc.yv, max_its_v[i], dims_v[i]);
        solve(A, b.zv, vc.zv, max_its_v[i], dims_v[i]);
        bool accurate = enforce_velocity_BC(BC, vc, accuracy_percent);
        if (accurate) {
            v_n += vc;
            const bool accurate2 =  enforce_velocity_BC(BC, v_n, accuracy_percent);
            if (accurate2) {
                std::cerr << "velocity solved using solver" << i << "\n";
                return;
            } else {
                v_n -= vc;
            }
        }
    }

    v_n += vc;
    std::cerr << "velocity is still not accurate even with most accurate solver\n";


    /*solve<normal_solve.max_iterations, normal_solve.dim>(A, b.xv, vc.xv);
    solve<normal_solve.max_iterations, normal_solve.dim>(A, b.yv, vc.yv);
    solve<normal_solve.max_iterations, normal_solve.dim>(A, b.zv, vc.zv);
    bool accurate = enforce_velocity_BC(BC, vc, accuracy_percent);
    if (accurate) {
        v_n += vc;
        const bool accurate2 = enforce_velocity_BC(BC, v_n, accuracy_percent);
        if (accurate2) {
            return;
        } else {
            v_n -= vc;
        }
    }

    solve<detailed_solve.max_iterations, detailed_solve.dim>(A, b.xv, vc.xv);
    solve<detailed_solve.max_iterations, detailed_solve.dim>(A, b.yv, vc.yv);
    solve<detailed_solve.max_iterations, detailed_solve.dim>(A, b.zv, vc.zv);
    accurate = enforce_velocity_BC(BC, vc, accuracy_percent);
    if (accurate) {
        v_n += vc;
        const bool accurate2 = enforce_velocity_BC(BC, v_n, accuracy_percent);
        if (accurate2) {
            return;
        } else {
            v_n -= vc;
        }
    }

    solve<very_detailed_solve.max_iterations, very_detailed_solve.dim>(A, b.xv, vc.xv);
    solve<very_detailed_solve.max_iterations, very_detailed_solve.dim>(A, b.yv, vc.yv);
    solve<very_detailed_solve.max_iterations, very_detailed_solve.dim>(A, b.zv, vc.zv);
    accurate = enforce_velocity_BC(BC, vc, accuracy_percent);
    if (accurate) {
        v_n += vc;
        const bool accurate2 = enforce_velocity_BC(BC, v_n, accuracy_percent);
        if (accurate2) {
            return;
        } else {
            std::cerr << "velocity is still not accurate even with most accurate solver\n";
        }
    }*/





}



#endif //CODE_SOLVER_HPP
