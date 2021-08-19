//
// Created by jacob on 12/8/21.
//

#ifndef CODE_BIG_MATRIX_HPP
#define CODE_BIG_MATRIX_HPP

#include <Eigen/Sparse>

template <unsigned N, unsigned M, unsigned P>
struct big_matrix{
    Eigen::SparseMatrix<double,Eigen::RowMajor> m;

    big_matrix() = delete;
    big_matrix(const unsigned no_data) {
        m.resize((N+1)*(M+1)*(P+1), (N+1)*(M+1)*(P+1));
        m.reserve(Eigen::VectorXi::Constant((N+1)*(M+1)*(P+1),no_data));
        //m.reserve( 8 * (N+1)*(M+1)*(P+1) ); //the amount of room to reserve --- overestimated
                                            // (7* is the default at non-boundary points
    };

    void add_elm(const unsigned i, const unsigned j, const unsigned k,
                 const unsigned x, const unsigned y, const unsigned z,
                 const double elm) {
        const auto row = get_index(i,j,k);
        const auto col = get_index(x,y,z);

        m.insert(row, col) = elm;
    }

private:
    [[nodiscard]] constexpr inline unsigned get_index(const unsigned i, const unsigned j, const unsigned k) const {
        return i + (N+1)*j + (N+1)*(M+1)*k;
    }

};


//for debugging
// - for the actual sparse matrix, just std::cout << mat.m
template <unsigned N, unsigned M, unsigned P>
std::ostream& operator << (std::ostream &out, const big_matrix<N,M,P> &mat) {
    Eigen::MatrixXd temp = mat.m;
    return out << temp;
}

#endif //CODE_BIG_MATRIX_HPP
