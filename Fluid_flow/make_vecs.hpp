//
// Created by jacob on 12/8/21.
//

#ifndef CODE_MAKE_VECS_HPP
#define CODE_MAKE_VECS_HPP

#include "../MyMath/calc.hpp"

//TODO : test to see if gives right output
template <unsigned N, unsigned M, unsigned P>
void make_b(big_vec<N,M,P,vec3> &b, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,vec3> &v_n1, const big_vec<N,M,P,double> &p, const big_vec<N,M,P,vec3> &bc) {
#pragma omp parallel for
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                const auto end_x = !v_n.has_left(i,j,k) || !v_n.has_right(i,j,k);
                const auto end_y = !v_n.has_down(i,j,k) || !v_n.has_up(i,j,k);
                const auto end_z = !v_n.has_front(i,j,k) || !v_n.has_back(i,j,k);
                if (end_x || end_y || end_z) {
                    b.add_elm(i,j,k,  bc(i,j,k));
                } else {
                    b.add_elm(i, j, k,
                              -3 / 2 * advection(v_n, i, j, k) + 1 / (2 * Re) * laplacian(v_n, i, j, k) +
                              v_n(i, j, k) / dt -
                              gradient(p, i, j, k) + 1 / 2 * advection(v_n1, i, j, k));
                }
            }
        }
    }


}

//TODO : test to see if gives right output
//for the first timestep
template <unsigned N, unsigned M, unsigned P>
void make_b_first(big_vec<N,M,P,vec3> &b, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,double> &p, const big_vec<N,M,P,vec3> &bc) {
#pragma omp parallel for
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                const auto end_x = !v_n.has_left(i,j,k) || !v_n.has_right(i,j,k);
                const auto end_y = !v_n.has_down(i,j,k) || !v_n.has_up(i,j,k);
                const auto end_z = !v_n.has_front(i,j,k) || !v_n.has_back(i,j,k);
                if (end_x || end_y || end_z) {
                    b.add_elm(i,j,k,  bc(i,j,k));
                } else {
                    b.add_elm(i, j, k,
                              -advection(v_n, i, j, k) + 1 / (2 * Re) * laplacian(v_n, i, j, k) + v_n(i, j, k) / dt -
                              gradient(p, i, j, k));
                }
            }
        }
    }
}

//TODO : test to see if gives right output
template <unsigned N, unsigned M, unsigned P>
void make_s(big_vec<N,M,P,double> &s, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,vec3> &v_n1, const big_vec<N,M,P,double> &p) {
#pragma omp parallel for
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                s(i,j,k) = divergence(v_n, i, j, k) / dt - 3/2 * divergence_advection(v_n, i,j,k) + 1/2 *
                        divergence_advection(v_n1, i,j,k) + 3/(2*Re) *
                        divergence_laplacian(v_n,i,j,k) - 1/(2*Re) * divergence_laplacian(v_n1, i,j,k) - laplacian(p, i,j,k);
            }
        }
    }
}


//TODO : test to see if gives right output
template <unsigned N, unsigned M, unsigned P>
void make_s_first(big_vec<N,M,P,double> &s, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,double> &p) {
#pragma omp parallel for
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                s(i,j,k) = divergence(v_n, i, j, k) / dt - divergence_advection(v_n,i,j,k) + 1/Re *
                                                                                                     divergence_laplacian(v_n,i,j,k)
                        -laplacian(p, i,j,k);
            }
        }
    }
}

#endif //CODE_MAKE_VECS_HPP
