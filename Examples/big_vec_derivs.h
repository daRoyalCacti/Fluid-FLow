//
// Created by jacob on 11/8/21.
//

#ifndef CODE_BIG_VEC_DERIVS_H
#define CODE_BIG_VEC_DERIVS_H

#include "../MyMath/big_vec.hpp"
#include "../MyMath/calc.hpp"
#include <iostream>

void big_veg_derivs_ex() {

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


    double x[nx][ny][nz], y[nx][ny][nz], z[nx][ny][nz];
    big_vec<nx-1, ny-1, nz-1, vec3> v(dx, dy, dz);

    for (unsigned i = 0; i < nx; i++) {
        for (unsigned j = 0; j < ny; j++) {
            for (unsigned k = 0; k < nz; k++) {
                x[i][j][k] = ex + dx*i;
                y[i][j][k] = ey + dy*j;
                z[i][j][k] = ez + dz*k;

                //v(i,j,k) = vec3(sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]));
                v.add_elm(i,j,k, sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]), sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]), sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]));
            }
        }
    }

    for (unsigned i = 0; i < nx; i++) {
        for (unsigned j = 0; j < ny; j++) {
            for (unsigned k = 0; k < nz; k++) {
                {   //d/dx
                    const auto true_deriv = cos(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<1, 0, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d/dx failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d/dy
                    const auto true_deriv = -sin(x[i][j][k])*sin(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 1, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d/dy failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //d/dz
                    const auto true_deriv = sin(x[i][j][k])*cos(y[i][j][k])*cos(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 0, 1>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d/dz failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }




                {   //d^2/dx^2
                    const auto true_deriv = -sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<2, 0, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^2/dx^2 failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^2/dy^2
                    const auto true_deriv = -sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 2, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^2/dy^2 failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //d^2/dz^2
                    const auto true_deriv = -sin(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 0, 2>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^2/dz^2 failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //d^2/dxdy
                    const auto true_deriv = -cos(x[i][j][k])*sin(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<1, 1, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^2/dxdy failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //d^2/dxdz
                    const auto true_deriv = cos(x[i][j][k])*cos(y[i][j][k])*cos(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<1, 0, 1>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^2/dxdz failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //d^2/dydz
                    const auto true_deriv = -sin(x[i][j][k])*sin(y[i][j][k])*cos(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 1, 1>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^2/dydz failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }




                {   //d^3/dx^3
                    const auto true_deriv = -cos(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<3, 0, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dx^3 failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dy^3
                    const auto true_deriv = sin(x[i][j][k])*sin(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 3, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dy^3 failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dz^3
                    const auto true_deriv = -sin(x[i][j][k])*cos(y[i][j][k])*cos(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 0, 3>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dz^3 failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }

                {   //d^3/dx^2dy
                    const auto true_deriv = sin(x[i][j][k])*sin(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<2, 1, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dx^2dy failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dy^2dx
                    const auto true_deriv = -cos(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<1, 2, 0>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dy^2dz failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dx^2dz
                    const auto true_deriv = -sin(x[i][j][k])*cos(y[i][j][k])*cos(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<2, 0, 1>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dx^2dz failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dz^2dx
                    const auto true_deriv = -cos(x[i][j][k])*cos(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<1, 0, 2>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dz^2dx failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dz^2dy
                    const auto true_deriv = sin(x[i][j][k])*sin(y[i][j][k])*sin(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 1, 2>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dz^2dy failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }
                {   //d^3/dy^2dz
                    const auto true_deriv = -sin(x[i][j][k])*cos(y[i][j][k])*cos(z[i][j][k]);
                    const auto calc_deriv = smart_deriv<0, 2, 1>(v, i, j, k).x();
                    if (abs(calc_deriv - true_deriv) > tol) {
                        std::cout << "d^3/dy^2dz failed at i = " << i << " j = " << j << " k = " << k << "\n";
                        std::cout << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
                    }
                }



            }
        }
    }

    std::cout << "checking complete\n";
}


#endif //CODE_BIG_VEC_DERIVS_H
