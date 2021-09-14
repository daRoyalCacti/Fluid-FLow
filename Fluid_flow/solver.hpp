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


//solves Ax=b for x
void solve(const big_matrix &A, const big_vec_d &b, big_vec_d &x) noexcept {
    viennacl::compressed_matrix<double>  vcl_sparse_matrix( b.g->size(), b.g->size() );
    //viennacl::coordinate_matrix<double>  vcl_sparse_matrix( (N+1)*(M+1)*(P+1), (N+1)*(M+1)*(P+1) );
    viennacl::vector<double> vcl_rhs( b.g->size() );

    //copying data into viennacl
    viennacl::copy(A.m, vcl_sparse_matrix);
    viennacl::copy(b.v, vcl_rhs);

    //setting up the solver
    //could probably mess with the parameters
    //viennacl::linalg::gmres_tag my_gmres_tag(1e-5, 100, 20); // up to 100 iterations, restart after 20 iterations
    //viennacl::linalg::gmres_tag my_gmres_tag(1e-10, 1000, 50); // up to 100 iterations, restart after 20 iterations
    viennacl::linalg::gmres_tag my_gmres_tag(1e-50, 1000, 20); // up to 100 iterations, restart after 20 iterations

    //solving
    viennacl::vector<double> res = viennacl::linalg::solve(vcl_sparse_matrix, vcl_rhs, my_gmres_tag);

    //copy result into x
    viennacl::copy(res, x.v);
}



#endif //CODE_SOLVER_HPP
