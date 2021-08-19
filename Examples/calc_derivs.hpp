//
// Created by jacob on 12/8/21.
//

#ifndef CODE_CALC_DERIVS_HPP
#define CODE_CALC_DERIVS_HPP

#include "../MyMath/calc.hpp"
#include "../MyMath/boundary.hpp"


namespace CD {
    template <unsigned N, unsigned M, unsigned P>
    void create_boundary_points(boundary_points<N,M,P> &bound) {
        for (unsigned i = 0; i <= N; i++) {
            for (unsigned k = 0; k <= P; k++) {
                bound(i,0,k).has_down = false;
                bound(i,M,k).has_up = false;
            }
        }
        for (unsigned i = 0; i <= N; i++) {
            for (unsigned j = 0; j <= M; j++) {
                bound(i,j,0).has_front = false;
                bound(i,j,P).has_back = false;
            }
        }
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                bound(0,j,k).has_left = false;
                bound(N,j,k).has_right = false;
            }
        }
    }
}

void calc_derivs_ex() {
    using namespace CD;

    constexpr unsigned nx = 10;
    constexpr unsigned ny = 16;
    constexpr unsigned nz = 12;

    constexpr double dx = 0.001;
    constexpr double dy = 0.0004;
    constexpr double dz = 0.0007;

    constexpr double ex = 3;
    constexpr double ey = 4;
    constexpr double ez = 5;

    constexpr double tol = 0.1;    //tolerance for if the derivative is correct

    std::cout << "testing to see if derivatives correct\n";

    boundary_points<nx-1,ny-1,nz-1> bound;
    create_boundary_points(bound);


    double xv[nx][ny][nz], yv[nx][ny][nz], zv[nx][ny][nz];
    big_vec<nx-1, ny-1, nz-1, vec3> v(dx, dy, dz, &bound);

    for (unsigned i = 0; i < nx; i++) {
        for (unsigned j = 0; j < ny; j++) {
            for (unsigned k = 0; k < nz; k++) {
                xv[i][j][k] = ex + dx*i;
                yv[i][j][k] = ey + dy*j;
                zv[i][j][k] = ez + dz*k;

                /*v(i,j,k) = vec3(cos(xv[i][j][k])*sin(yv[i][j][k])*sin(zv[i][j][k]),
                                sin(xv[i][j][k])*cos(yv[i][j][k])*sin(zv[i][j][k]),
                                sin(xv[i][j][k])*sin(yv[i][j][k])*cos(zv[i][j][k]));*/
                v.add_elm(i,j,k,
                          cos(xv[i][j][k])*sin(yv[i][j][k])*sin(zv[i][j][k]),
                          sin(xv[i][j][k])*cos(yv[i][j][k])*sin(zv[i][j][k]),
                          sin(xv[i][j][k])*sin(yv[i][j][k])*cos(zv[i][j][k]) );
            }
        }
    }

    for (unsigned i = 0; i < nx; i++) {
        for (unsigned j = 0; j < ny; j++) {
            for (unsigned k = 0; k < nz; k++) {
                {   //advection
                    const auto x = xv[i][j][k];
                    const auto y = yv[i][j][k];
                    const auto z = zv[i][j][k];

                    const auto true_deriv = vec3(cos(x) * pow(cos(z),2) *sin(x) *pow(sin(y),2) + cos(x) *pow(cos(y),2)* sin(x) *pow(sin(z),2) - cos(x) *sin(x)* pow(sin(y),2) *pow(sin(z),2),
                                                 cos(y)* pow(cos(z),2) *pow(sin(x),2)* sin(y) + pow(cos(x),2)* cos(y) *sin(y) *pow(sin(z),2) - cos(y) *pow(sin(x),2) *sin(y) *pow(sin(z),2),
                                                 pow(cos(y),2) *cos(z)* pow(sin(x),2) *sin(z) + pow(cos(x),2) *cos(z) *pow(sin(y),2) *sin(z) - cos(z) *pow(sin(x),2)* pow(sin(y),2) *sin(z));
                    const auto calc_deriv = advection(v, i, j, k);
                    if ((calc_deriv - true_deriv).length_squared() > tol) {
                        std::cerr << "advection failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //lapacian
                    const auto x = xv[i][j][k];
                    const auto y = yv[i][j][k];
                    const auto z = zv[i][j][k];

                    const auto true_deriv = vec3(-3*cos(x)*sin(y)*sin(z), -3*cos(y)*sin(x)*sin(z), -3*cos(z)*sin(x)*sin(y));
                    const auto calc_deriv = laplacian(v, i, j, k);
                    if ((calc_deriv - true_deriv).length_squared() > tol) {
                        std::cerr << "lapacian failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //Divergence
                    const auto x = xv[i][j][k];
                    const auto y = yv[i][j][k];
                    const auto z = zv[i][j][k];

                    const auto true_deriv =-3*sin(x)*sin(y)*sin(z);
                    const auto calc_deriv = divergence(v, i, j, k);
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cerr << "Divergence failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //Divergence of advection
                    const auto x = xv[i][j][k];
                    const auto y = yv[i][j][k];
                    const auto z = zv[i][j][k];

                    const auto true_deriv =2* pow(cos(y),2)* pow(cos(z),2) *pow(sin(x),2) + 2 *pow(cos(x),2) *pow(cos(z),2)* pow(sin(y),2) -
                                           3* pow(cos(z),2) *pow(sin(x),2) *pow(sin(y),2) + 2 *pow(cos(x),2) *pow(cos(y),2) *pow(sin(z),2) -
                                           3 *pow(cos(y),2)* pow(sin(x),2)* pow(sin(z),2) - 3 *pow(cos(x),2) *pow(sin(y),2) *pow(sin(z),2) +
                                           pow(sin(x),2)* pow(sin(y),2) *pow(sin(z),2);
                    const auto calc_deriv = divergence_advection(v, i, j, k);
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cerr << "Divergence of advection failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //Divergence of laplacian
                    const auto x = xv[i][j][k];
                    const auto y = yv[i][j][k];
                    const auto z = zv[i][j][k];

                    const auto true_deriv =9 * sin(x) * sin(y) * sin(z);
                    const auto calc_deriv = divergence_laplacian(v, i, j, k);
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cerr << "Divergence of laplacian failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

            }
        }
    }


    std::cout << "testing complete\n";
}

#endif //CODE_CALC_DERIVS_HPP
