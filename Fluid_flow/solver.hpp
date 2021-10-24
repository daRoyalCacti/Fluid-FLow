//
// Created by jacob on 15/8/21.
//

#ifndef CODE_SOLVER_HPP
#define CODE_SOLVER_HPP

#define CL_TARGET_OPENCL_VERSION 300    //the lastest version supported atm
#define VIENNACL_WITH_EIGEN 1

//#define BICGSTAB

#define IHU0
//#define JACOBI

#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <eigen3/unsupported/Eigen/IterativeSolvers>

#ifndef BICGSTAB
#include <viennacl/linalg/gmres.hpp>
#else
#include <viennacl/linalg/bicgstab.hpp>
#endif

#ifdef IHU0
#include <viennacl/linalg/detail/ilu/ilu0.hpp>
#endif

//#include <viennacl/linalg/ichol.hpp>
//#include <viennacl/linalg/detail/ilu/block_ilu.hpp>
#ifdef JACOBI
#include <viennacl/linalg/jacobi_precond.hpp>
#endif

#include "../MyMath/big_matrix.hpp"
#include "../MyMath/big_vec.hpp"
#include "update_vecs.hpp"


/*
constexpr std::array<unsigned, 5> max_its_v = {300, 500, 1000, 5000, 10000,};
constexpr std::array<unsigned, 5>  dims_v = {20, 20, 20, 20, 25};
constexpr unsigned no_solver_choices_v = max_its_v.size();*/

//can solve v with a tolerance of 10000, with 1 iteration and 1 dimension at a tolerance of 0.01
// - was no faster than with normal variables --- think all the time is spent uploading to device
constexpr std::array<double, 5> tols_v = {1e-3,1e-4, 1e-5, 1e-6, 1e-7};
constexpr std::array<unsigned, 5> max_its_v = {300, 500, 1000, 5000, 10000};
constexpr unsigned no_solver_choices_v = max_its_v.size();

constexpr unsigned krylov_dim = 10;
constexpr std::array<double, 4> tols_p = { 5e-10, 2e-10, 1e-10, 1e-20};
//constexpr std::array<unsigned, 4> max_its_p = {5000, 5000, 5000, 10000};
constexpr std::array<unsigned, 4> max_its_p = {50000, 50000, 50000, 100000};
constexpr std::array<unsigned, 4> dims_p    = {krylov_dim, krylov_dim, krylov_dim, krylov_dim};
constexpr unsigned no_solver_choices_p = max_its_p.size();

/*constexpr std::array<double, 10> tols_p = { 1e-1, 1e-2, 1e-3,1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
constexpr std::array<unsigned, 10> max_its_p = {5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000};
constexpr std::array<unsigned, 10> dims_p    = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
constexpr unsigned no_solver_choices_p = max_its_p.size();*/

/*
constexpr std::array<unsigned, 3> max_its_p = {3000, 5000, 10000};
constexpr std::array<unsigned, 3> dims_p    = {20, 20, 25};
constexpr unsigned no_solver_choices_p = max_its_p.size();*/


//solves Ax=b for x
//template <unsigned it, unsigned dim>
void solve_GPU(const big_matrix &A, const big_vec_d &b, big_vec_d &x, const unsigned it, const unsigned dim, const double tol = 1e-100) noexcept {
    viennacl::compressed_matrix<double>  vcl_sparse_matrix( b.size(), b.size() );
    viennacl::vector<double> vcl_rhs( b.size() );

    //copying data into viennacl
    viennacl::copy(A.m, vcl_sparse_matrix);
    viennacl::copy(b.v, vcl_rhs);

    //setting up the solver
    //could probably mess with the parameters
#ifdef VIENNACL_WITH_OPENCL
#ifndef BICGSTAB

viennacl::linalg::gmres_tag my_gmres_tag(tol, it, dim);

#ifdef IHU0
viennacl::linalg::ilu0_tag ilu0_config;
ilu0_config.use_level_scheduling(true);
viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<double> > vcl_ilut(vcl_sparse_matrix, ilu0_config);
viennacl::vector<double> res = viennacl::linalg::solve(vcl_sparse_matrix, vcl_rhs, my_gmres_tag, vcl_ilut);
#else
#ifdef JACOBI
viennacl::linalg::jacobi_precond< viennacl::compressed_matrix<double> > vcl_jacobi(vcl_sparse_matrix, viennacl::linalg::jacobi_tag());
viennacl::vector<double> res = viennacl::linalg::solve(vcl_sparse_matrix, vcl_rhs, my_gmres_tag, vcl_jacobi);
#else
viennacl::vector<double> res = viennacl::linalg::solve(vcl_sparse_matrix, vcl_rhs, my_gmres_tag);
#endif
#endif
//viennacl::linalg::block_ilu_precond< viennacl::compressed_matrix<double>, viennacl::linalg::ilu0_tag> vcl_block_ilu0(vcl_sparse_matrix, viennacl::linalg::ilu0_tag());
//viennacl::linalg::row_scaling< viennacl::compressed_matrix<double> > vcl_row_scaling(vcl_sparse_matrix, viennacl::linalg::row_scaling_tag());





/*viennacl::vector<double> x_guess( x.size() );
viennacl::copy(x.v, x_guess);
viennacl::linalg::gmres_tag my_gmres_tag(1e-100, it, dim);
viennacl::linalg::gmres_solver<viennacl::vector<double>> s(my_gmres_tag);
s.set_initial_guess(x_guess);
viennacl::vector<double> res = s(vcl_sparse_matrix, vcl_rhs);*/
#else
viennacl::vector<double> x_guess( x.size() );
viennacl::copy(x.v, x_guess);
viennacl::linalg::bicgstab_tag my_bicgstab_tag(1e-100, it, dim);
viennacl::linalg::bicgstab_solver<viennacl::vector<double>> s(my_bicgstab_tag);
s.set_initial_guess(x_guess);
viennacl::vector<double> res = s(vcl_sparse_matrix, vcl_rhs);
#endif
#endif
/*#ifdef VIENNACL_WITH_OPENMP
viennacl::linalg::gmres_tag my_gmres_tag(1e-5, it, dim);
#endif*/


//copy result into x
viennacl::copy(res, x.v);
}

void solveCPU(const big_matrix &A, const big_vec_d &b, big_vec_d &x, const unsigned it, const double tol = 1e-100) noexcept {
    Eigen::GMRES<Eigen::SparseMatrix<double,Eigen::RowMajor>> solver(A.m);
    solver.setTolerance(tol);
    solver.setMaxIterations(it);
    x.v = solver.solveWithGuess(b.v, x.v);
    //x.v = solver.solve(b.v);
}


void solve_pressure(const boundary_conditions& BC, const big_matrix &Q, const big_vec_d &s, big_vec_d &p_c, big_vec_d &p, const double accuracy_percent) {
    for (unsigned i; i < no_solver_choices_p; i++) {
        auto p_cpy = p;
        solve_GPU(Q, s, p_c, max_its_p[i], dims_p[i], tols_p[i]);
        bool accurate = update_pressure_BC<true>(BC, p_c, accuracy_percent);
        if (accurate) {
            p_cpy += p_c;
            const bool accurate2 = update_pressure_BC<true>(BC, p_cpy, accuracy_percent);
            if (accurate2) {
                p = p_cpy;
                std::cerr << "pressure solved using solver " << i << "\n";
                std::cerr << "\tits:" << max_its_p[i] << " dim:" << dims_p[i] << " tol:" << tols_p[i] << "\n";
                return;
            } else {
                std::cerr << "solver " << i << " failed\n";
                std::cerr << "\tits:" << max_its_p[i] << " dim:" << dims_p[i] << " tol:" << tols_p[i] << "\n";
            }
        } else {
            std::cerr << "solver " << i << " failed\n";
            std::cerr << "\tits:" << max_its_p[i] << " dim:" << dims_p[i] << " tol:" << tols_p[i] << "\n";
        }
    }

    std::cerr << "pressure is still not accurate even with most accurate solver\n";
    const auto i = no_solver_choices_p-1;
    solve_GPU(Q, s, p_c, max_its_p[i], dims_p[i]);
    update_pressure_BC<false>(BC, p_c, accuracy_percent);
    p += p_c;
    update_pressure_BC<false>(BC, p, accuracy_percent);

}

void solve_velocity(const boundary_conditions& BC, const big_matrix &A, const big_vec_v &b, big_vec_v &vc, big_vec_v &v_n, const double accuracy_percent) {
    for (unsigned i; i < no_solver_choices_v; i++) {
        auto v_n_cpy = v_n;
        solveCPU(A, b.xv, vc.xv, max_its_v[i], tols_v[i]);
        solveCPU(A, b.yv, vc.yv, max_its_v[i], tols_v[i]);
        solveCPU(A, b.zv, vc.zv, max_its_v[i], tols_v[i]);
        bool accurate = enforce_velocity_correction_BC<true>(BC, vc, accuracy_percent);
        if (accurate) {
            v_n_cpy += vc;
            const bool accurate2 =  enforce_velocity_BC<true>(BC, v_n_cpy, accuracy_percent);
            if (accurate2) {
                v_n = v_n_cpy;
                std::cerr << "velocity solved using solver " << i << "\n";
                std::cerr << "\tits:" << max_its_v[i]  << " tol:" << tols_v[i] << "\n";
                return;
            } else {
                std::cerr << "solver " << i << " failed\n";
                std::cerr << "\tits:" << max_its_v[i] << " tol:" << tols_v[i] << "\n";
            }
        } else {
            std::cerr << "solver " << i << " failed\n";
            std::cerr << "\tits:" << max_its_v[i]  << " tol:" << tols_v[i] << "\n";
        }
    }


    std::cerr << "velocity is still not accurate even with most accurate solver\n";
    const auto i = no_solver_choices_v-1;

    solveCPU(A, b.xv, vc.xv, max_its_v[i]);
    solveCPU(A, b.yv, vc.yv, max_its_v[i]);
    solveCPU(A, b.zv, vc.zv, max_its_v[i]);
    enforce_velocity_correction_BC<false>(BC, vc, accuracy_percent);
    v_n += vc;
    enforce_velocity_BC<false>(BC, v_n, accuracy_percent);


}



#endif //CODE_SOLVER_HPP
