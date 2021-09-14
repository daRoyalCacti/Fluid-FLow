//
// Created by jacob on 12/8/21.
//

#ifndef CODE_CALC_DERIVS_HPP
#define CODE_CALC_DERIVS_HPP

#include "../MyMath/calc.hpp"
#include "../MyMath/boundary.hpp"


void calc_derivs_ex() {

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

    boundary_conditions BC( dx, dy, dz, nx*dx, ny*dy, nz*dz, ex, ey, ez);
    big_vec_v v(BC);

    for (unsigned i = 0; i < v.size(); i++) {
        const auto x = BC.global_grid.x[i];
        const auto y = BC.global_grid.y[i];
        const auto z = BC.global_grid.z[i];
        v.add_elm(i, cos(x)*sin(y)*sin(z),sin(x)*cos(y)*sin(z),sin(x)*sin(y)*cos(z) );

    }

    for (unsigned i = 0; i < v.size(); i++) {
        const auto x = BC.global_grid.x[i];
        const auto y = BC.global_grid.y[i];
        const auto z = BC.global_grid.z[i];
        {   //advection
            const auto true_deriv = vec3(cos(x) * pow(cos(z),2) *sin(x) *pow(sin(y),2) + cos(x) *pow(cos(y),2)* sin(x) *pow(sin(z),2) - cos(x) *sin(x)* pow(sin(y),2) *pow(sin(z),2),
                                         cos(y)* pow(cos(z),2) *pow(sin(x),2)* sin(y) + pow(cos(x),2)* cos(y) *sin(y) *pow(sin(z),2) - cos(y) *pow(sin(x),2) *sin(y) *pow(sin(z),2),
                                         pow(cos(y),2) *cos(z)* pow(sin(x),2) *sin(z) + pow(cos(x),2) *cos(z) *pow(sin(y),2) *sin(z) - cos(z) *pow(sin(x),2)* pow(sin(y),2) *sin(z));
            const auto calc_deriv = advection(v, i);
            if ((calc_deriv - true_deriv).length_squared() > tol) {
                std::cerr << "advection failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //lapacian
            const auto true_deriv = vec3(-3*cos(x)*sin(y)*sin(z), -3*cos(y)*sin(x)*sin(z), -3*cos(z)*sin(x)*sin(y));
            const auto calc_deriv = laplacian(v, i);
            if ((calc_deriv - true_deriv).length_squared() > tol) {
                std::cerr << "lapacian failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //Divergence
            const auto true_deriv =-3*sin(x)*sin(y)*sin(z);
            const auto calc_deriv = divergence(v, i);
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "Divergence failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //Divergence of advection
            const auto true_deriv =2* pow(cos(y),2)* pow(cos(z),2) *pow(sin(x),2) + 2 *pow(cos(x),2) *pow(cos(z),2)* pow(sin(y),2) -
                    3* pow(cos(z),2) *pow(sin(x),2) *pow(sin(y),2) + 2 *pow(cos(x),2) *pow(cos(y),2) *pow(sin(z),2) -
                    3 *pow(cos(y),2)* pow(sin(x),2)* pow(sin(z),2) - 3 *pow(cos(x),2) *pow(sin(y),2) *pow(sin(z),2) +
                    pow(sin(x),2)* pow(sin(y),2) *pow(sin(z),2);
            const auto calc_deriv = divergence_advection(v, i);
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "Divergence of advection failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //Divergence of laplacian
            const auto true_deriv = 9 * sin(x) * sin(y) * sin(z);
            const auto calc_deriv = divergence_laplacian(v, i);
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "Divergence of laplacian failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

    }


    std::cout << "testing complete\n";
}

#endif //CODE_CALC_DERIVS_HPP
