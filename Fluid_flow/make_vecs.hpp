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
                if (v_n.is_boundary(i,j,k)) {
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
                if (v_n.is_boundary(i,j,k)) {
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
#ifdef pressure_BC
void make_s(big_vec<N,M,P,double> &s, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n,
        const big_vec<N,M,P,vec3> &v_n1, const big_vec<N,M,P,double> &p, const boundary_normals<N,M,P> &norms) {
#else
void make_s(big_vec<N,M,P,double> &s, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,vec3> &v_n1, const big_vec<N,M,P,double> &p) {
#endif
#pragma omp parallel for
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
#ifdef pressure_BC
                if (p.is_boundary(i,j,k)) {
                    const auto norm = norms.normal(i,j,k);
                    const auto nx = norm.x();
                    const auto ny = norm.y();
                    const auto nz = norm.z();

                    const auto dx = p.dx;
                    const auto dy = p.dy;
                    const auto dz = p.dz;

                    //picking the direction
                    unsigned big_dir = 0;
                    if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                        big_dir = 1;
                    }
                    if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                        big_dir = 2;
                    }


                    if (big_dir == 0) { //x direction biggest
                        if (!p.has_left(i,j,k)) {   //forward difference
                            s(i,j,k) = ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i+2,j,k)/3 + 4*p(i+1,j,k)/3;
                        } else if (!p.has_right(i,j,k)) {   //backward difference
                            s(i,j,k) =-ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i-2,j,k)/3 + 4*p(i-1,j,k)/3;
                        } else {    //central difference
                            s(i,j,k) = -ny/nx * smart_deriv<0,1,0>(p, i,j,k)*dx - nz/nx * smart_deriv<0,0,1>(p, i,j,k)*dx + p(i-1,j,k);
                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p.has_down(i, j, k)) { //forward difference
                            s(i,j,k) = nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j+2,k)/3 + 4*p(i,j+1,k)/3;
                        } else if (!p.has_up(i, j, k)) {  //backward difference
                            s(i,j,k) =-nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j-2,k)/3 + 4*p(i,j-1,k)/3;
                        } else {  //central difference
                            s(i,j,k) = -nx/ny * smart_deriv<1,0,0>(p, i,j,k)*dy - nz/ny * smart_deriv<0,0,1>(p, i,j,k)*dy + p(i,j-1,k);
                        }
                    } else {    //z direction
#ifndef NDEBUG
                        if (big_dir > 2) {
                            std::cerr << "the biggest direction cannot be larger than 2\n";
                        }
#endif
                        if (!p.has_front(i,j,k)) {  //forward difference
                            s(i,j,k) = ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k+2)/3 + 4*p(i,j,k+1)/3;
                        } else if (!p.has_back(i,j,k)) { //backward difference
                            s(i,j,k) =-ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k-2)/3 + 4*p(i,j,k-1)/3;
                        } else {
                            s(i,j,k) = -ny/nz * smart_deriv<0,1,0>(p, i,j,k)*dz - nx/nz * smart_deriv<1,0,0>(p, i,j,k)*dz + p(i,j,k-1);
                        }
                    }

                } else {
#endif

                s(i,j,k) = divergence(v_n, i, j, k) / dt - 3/2 * divergence_advection(v_n, i,j,k) + 1/2 *
                        divergence_advection(v_n1, i,j,k) + 3/(2*Re) *
                        divergence_laplacian(v_n,i,j,k) - 1/(2*Re) * divergence_laplacian(v_n1, i,j,k) - laplacian(p, i,j,k);
#ifdef pressure_BC
                }
#endif
            }
        }
    }
}


//TODO : test to see if gives right output
template <unsigned N, unsigned M, unsigned P>
#ifdef pressure_BC
void make_s_first(big_vec<N,M,P,double> &s, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,double> &p, const boundary_normals<N,M,P> &norms) {
#else
void make_s_first(big_vec<N,M,P,double> &s, const double Re, const double dt, const big_vec<N,M,P,vec3> &v_n, const big_vec<N,M,P,double> &p) {
#endif
#pragma omp parallel for
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
#ifdef pressure_BC
                if (p.is_boundary(i,j,k)) {

                    const auto norm = norms.normal(i,j,k);
                    const auto nx = norm.x();
                    const auto ny = norm.y();
                    const auto nz = norm.z();

                    const auto dx = p.dx;
                    const auto dy = p.dy;
                    const auto dz = p.dz;

                    //picking the direction
                    unsigned big_dir = 0;
                    if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                        big_dir = 1;
                    }
                    if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                        big_dir = 2;
                    }


                    if (big_dir == 0) { //x direction biggest
                        if (!p.has_left(i,j,k)) {   //forward difference
                            s(i,j,k) = ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i+2,j,k)/3 + 4*p(i+1,j,k)/3;
                        } else if (!p.has_right(i,j,k)) {   //backward difference
                            s(i,j,k) =-ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i-2,j,k)/3 + 4*p(i-1,j,k)/3;
                        } else {    //central difference
                            s(i,j,k) = -ny/nx * smart_deriv<0,1,0>(p, i,j,k)*dx - nz/nx * smart_deriv<0,0,1>(p, i,j,k)*dx + p(i-1,j,k);
                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p.has_down(i, j, k)) { //forward difference
                            s(i,j,k) = nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j+2,k)/3 + 4*p(i,j+1,k)/3;
                        } else if (!p.has_up(i, j, k)) {  //backward difference
                            s(i,j,k) =-nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j-2,k)/3 + 4*p(i,j-1,k)/3;
                        } else {  //central difference
                            s(i,j,k) = -nx/ny * smart_deriv<1,0,0>(p, i,j,k)*dy - nz/ny * smart_deriv<0,0,1>(p, i,j,k)*dy + p(i,j-1,k);
                        }
                    } else {    //z direction
#ifndef NDEBUG
                        if (big_dir > 2) {
                            std::cerr << "the biggest direction cannot be larger than 2\n";
                        }
#endif
                        if (!p.has_front(i,j,k)) {  //forward difference
                            s(i,j,k) = ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k+2)/3 + 4*p(i,j,k+1)/3;
                        } else if (!p.has_back(i,j,k)) { //backward difference
                            s(i,j,k) =-ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k-2)/3 + 4*p(i,j,k-1)/3;
                        } else {
                            s(i,j,k) = -ny/nz * smart_deriv<0,1,0>(p, i,j,k)*dz - nx/nz * smart_deriv<1,0,0>(p, i,j,k)*dz + p(i,j,k-1);
                        }
                    }

                } else {
#endif
                    s(i, j, k) = divergence(v_n, i, j, k) / dt - divergence_advection(v_n, i, j, k) + 1 / Re *
                                  divergence_laplacian(v_n, i, j,k) - laplacian(p, i, j, k);
#ifdef pressure_BC
                }
#endif
            }
        }
    }
}

#endif //CODE_MAKE_VECS_HPP
