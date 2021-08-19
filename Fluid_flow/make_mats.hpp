//
// Created by jacob on 13/8/21.
//

#ifndef CODE_MAKE_MATS_HPP
#define CODE_MAKE_MATS_HPP

#include "../MyMath/big_matrix.hpp"
#include "../MyMath/big_vec.hpp"

/*
template <unsigned N, unsigned M, unsigned P>
void make_Q(big_matrix<N,M,P> &Q, const big_vec<N,M,P,double> &p) {
    const auto dxdx = p.dx*p.dx;
    const auto dydy = p.dy*p.dy;
    const auto dzdz = p.dz*p.dz;

    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                double diag = 0.0;

                //x axis checks
                if (!p.has_left(i,j,k)) {
                    diag += 2/dxdx;

                    Q.add_elm(i,j,k,  i+1,j,k,  1/dxdx);
                } else if (!p.has_right(i,j,k)) {
                    diag += 2/dxdx;

                    Q.add_elm(i,j,k,  i-1,j,k,  1/dxdx);
                } else {
                    diag += -2/dxdx;

                    Q.add_elm(i,j,k,  i+1,j,k,  1/dxdx);
                    Q.add_elm(i,j,k,  i-1,j,k,  1/dxdx);
                }

                //y axis checks
                if (!p.has_down(i,j,k)) {
                    diag += 2/dydy;

                    Q.add_elm(i,j,k,  i,j+1,k,  1/dydy);
                } else if (!p.has_up(i,j,k)) {
                    diag += 2/dydy;

                    Q.add_elm(i,j,k,  i,j-1,k,  1/dydy);
                } else {
                    diag += -2/dydy;

                    Q.add_elm(i,j,k,  i,j+1,k,  1/dydy);
                    Q.add_elm(i,j,k,  i,j-1,k,  1/dydy);
                }

                //z axis checks
                if (!p.has_front(i,j,k)) {
                    diag += 2/dzdz;

                    Q.add_elm(i,j,k,  i,j,k+1,  1/dzdz);
                } else if (!p.has_back(i,j,k)) {
                    diag += 2/dzdz;

                    Q.add_elm(i,j,k,  i,j,k-1,  1/dzdz);
                } else {
                    diag += -2/dzdz;

                    Q.add_elm(i,j,k,  i,j,k+1,  1/dzdz);
                    Q.add_elm(i,j,k,  i,j,k-1,  1/dzdz);
                }

                Q.add_elm(i,j,k,  i,j,k,  diag);
            }
        }
    }

}
 */

//p to store dx,dy,dz and BC
template <unsigned N, unsigned M, unsigned P>
void make_Q(big_matrix<N,M,P> &Q, const big_vec<N,M,P,double> &p) {
    const auto dxdx = p.dx*p.dx;
    const auto dydy = p.dy*p.dy;
    const auto dzdz = p.dz*p.dz;

//#pragma omp parallel for
    //shared(Q, p, dxdx, dydy, dzdz) default(none)
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                double diag = 0.0;

                //x axis checks
                if (!p.has_left(i,j,k)) {
                    diag += 2/dxdx;

                    Q.add_elm(i,j,k,  i+3,j,k, -1/dxdx);
                    Q.add_elm(i,j,k,  i+2,j,k, 4/dxdx);
                    Q.add_elm(i,j,k,  i+1,j,k, -5/dxdx);
                } else if (!p.has_right(i,j,k)) {
                    diag += 2/dxdx;

                    Q.add_elm(i,j,k,  i-1,j,k, -5/dxdx);
                    Q.add_elm(i,j,k,  i-2,j,k, 4/dxdx);
                    Q.add_elm(i,j,k,  i-3,j,k, -1/dxdx);
                } else {
                    diag += -2/dxdx;

                    Q.add_elm(i,j,k,  i+1,j,k,  1/dxdx);
                    Q.add_elm(i,j,k,  i-1,j,k,  1/dxdx);
                }

                //y axis checks
                if (!p.has_down(i,j,k)) {
                    diag += 2/dydy;

                    Q.add_elm(i,j,k,  i,j+3,k, -1/dydy);
                    Q.add_elm(i,j,k,  i,j+2,k, 4/dydy);
                    Q.add_elm(i,j,k,  i,j+1,k, -5/dydy);
                } else if (!p.has_up(i,j,k)) {
                    diag += 2/dydy;

                    Q.add_elm(i,j,k,  i,j-1,k, -5/dydy);
                    Q.add_elm(i,j,k,  i,j-2,k, 4/dydy);
                    Q.add_elm(i,j,k,  i,j-3,k, -1/dydy);
                } else {
                    diag += -2/dydy;

                    Q.add_elm(i,j,k,  i,j+1,k,  1/dydy);
                    Q.add_elm(i,j,k,  i,j-1,k,  1/dydy);
                }

                //z axis checks
                if (!p.has_front(i,j,k)) {
                    diag += 2/dzdz;

                    Q.add_elm(i,j,k,  i,j,k+3, -1/dzdz);
                    Q.add_elm(i,j,k,  i,j,k+2, 4/dzdz);
                    Q.add_elm(i,j,k,  i,j,k+1, -5/dzdz);
                } else if (!p.has_back(i,j,k)) {
                    diag += 2/dzdz;

                    Q.add_elm(i,j,k,  i,j,k-1, -5/dzdz);
                    Q.add_elm(i,j,k,  i,j,k-2, 4/dzdz);
                    Q.add_elm(i,j,k,  i,j,k-3, -1/dzdz);
                } else {
                    diag += -2/dzdz;

                    Q.add_elm(i,j,k,  i,j,k+1,  1/dzdz);
                    Q.add_elm(i,j,k,  i,j,k-1,  1/dzdz);
                }

                Q.add_elm(i,j,k,  i,j,k,  diag);

            }
        }
    }
}



template <unsigned N, unsigned M, unsigned P>
void make_A(big_matrix<N,M,P> &A,  const big_vec<N,M,P,vec3> &v, const double dt, const double Re) {
    const auto Rdxdx = Re*v.dx*v.dx;
    const auto Rdydy = Re*v.dy*v.dy;
    const auto Rdzdz = Re*v.dz*v.dz;

//#pragma omp parallel for
    //shared(A, v, dt, Re, Rdxdx, Rdydy, Rdzdz) default(none)
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                if (v.is_boundary(i,j,k)) {
                    A.add_elm(i,j,k,  i,j,k,  1);
                } else {
                    A.add_elm(i,j,k,  i,j,k,  1/dt + 1/Rdxdx + 1/Rdydy + 1/Rdzdz );
                    A.add_elm(i,j,k,  i+1,j,k,  -1/(2*Rdxdx));
                    A.add_elm(i,j,k,  i-1,j,k,  -1/(2*Rdxdx));
                    A.add_elm(i,j,k,  i,j+1,k,  -1/(2*Rdydy));
                    A.add_elm(i,j,k,  i,j-1,k,  -1/(2*Rdydy));
                    A.add_elm(i,j,k,  i,j,k+1,  -1/(2*Rdzdz));
                    A.add_elm(i,j,k,  i,j,k-1,  -1/(2*Rdzdz));
                }
            }
        }
    }
}

#endif //CODE_MAKE_MATS_HPP
