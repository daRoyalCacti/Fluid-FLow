//
// Created by jacob on 12/8/21.
//

#ifndef CODE_BIG_MATRIX_HPP
#define CODE_BIG_MATRIX_HPP

#include <Eigen/Sparse>
#include <iostream>

struct big_matrix{
    Eigen::SparseMatrix<double,Eigen::RowMajor> m;

    big_matrix() = delete;
    big_matrix(const boundary_conditions &b, const unsigned no_data) noexcept {
        m.resize(static_cast<long>(b.global_grid.size()), static_cast<long>(b.global_grid.size())  );
        m.reserve(Eigen::VectorXi::Constant(static_cast<long>(b.global_grid.size() ), static_cast<int>(no_data) ) );
        //m.reserve( 8 * (N+1)*(M+1)*(P+1) ); //the amount of room to reserve --- overestimated
                                            // (7* is the default at non-boundary points
    };

    void add_elm(const unsigned index1, const unsigned index2, const double elm) noexcept {
        m.insert(index1, index2) = elm;
    }


};


//for debugging
// - for the actual sparse matrix, just std::cout << mat.m
std::ostream& operator << (std::ostream &out, const big_matrix &mat) {
    Eigen::MatrixXd temp = mat.m;
    return out << temp;
}

#endif //CODE_BIG_MATRIX_HPP
