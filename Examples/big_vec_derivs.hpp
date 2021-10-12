//
// Created by jacob on 11/8/21.
//

#ifndef CODE_BIG_VEC_DERIVS_HPP
#define CODE_BIG_VEC_DERIVS_HPP

#include "../MyMath/big_vec.hpp"
#include "../MyMath/calc.hpp"
#include "../MyMath/boundary.hpp"
#include <iostream>
#include "../Fluid_flow/boundary_conditions.hpp"

//#define BVED_DLOG   //bvec example derivating detailed logging

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

    boundary_conditions BC(dx, dy, dz, nx*dx, ny*dy, nz*dz, ex, ey, ez);
    big_vec_v v(BC);


    for (unsigned i = 0; i < v.size(); i++) {
        //std::cerr << i << "/" << v.size() << "\n";
        const auto x = BC.global_grid.x[i];
        const auto y = BC.global_grid.y[i];
        const auto z = BC.global_grid.z[i];
        v.add_elm(i, sin(x)*cos(y)*sin(z), sin(x)*cos(y)*sin(z), sin(x)*cos(y)*sin(z) );
    }

    for (unsigned i = 0; i < v.size(); i++) {
#ifdef BVED_DLOG
        std::cerr << i << "/" << v.size() << "\n";
#endif
        const auto x = BC.global_grid.x[i];
        const auto y = BC.global_grid.y[i];
        const auto z = BC.global_grid.z[i];
        {   //d/dx
#ifdef BVED_DLOG
            std::cerr << "\td/dx\n";
            #endif
            const auto true_deriv = cos(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<1, 0, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d/dx failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d/dy
            #ifdef BVED_DLOG
            std::cerr << "\td/dy\n";
            #endif
            const auto true_deriv = -sin(x)*sin(y)*sin(z);
            const auto calc_deriv = smart_deriv<0, 1, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d/dy failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //d/dz
            #ifdef BVED_DLOG
            std::cerr << "\td/dz\n";
            #endif
            const auto true_deriv = sin(x)*cos(y)*cos(z);
            const auto calc_deriv = smart_deriv<0, 0, 1>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d/dz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }




        {   //d^2/dx^2
            #ifdef BVED_DLOG
            std::cerr << "\td^2/dx^2\n";
            #endif
            const auto true_deriv = -sin(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<2, 0, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^2/dx^2 failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^2/dy^2
            #ifdef BVED_DLOG
            std::cerr << "\td^2/dx^2\n";
            #endif
            const auto true_deriv = -sin(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<0, 2, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^2/dy^2 failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //d^2/dz^2
            #ifdef BVED_DLOG
            std::cerr << "\td^2/dz^2\n";
            #endif
            const auto true_deriv = -sin(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<0, 0, 2>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^2/dz^2 failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //d^2/dxdy
            #ifdef BVED_DLOG
            std::cerr << "\td^2/dxdy\n";
            #endif
            const auto true_deriv = -cos(x)*sin(y)*sin(z);
            const auto calc_deriv = smart_deriv<1, 1, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^2/dxdy failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //d^2/dxdz
            #ifdef BVED_DLOG
            std::cerr << "\td^2/dxdz\n";
            #endif
            const auto true_deriv = cos(x)*cos(y)*cos(z);
            const auto calc_deriv = smart_deriv<1, 0, 1>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^2/dxdz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //d^2/dydz
            #ifdef BVED_DLOG
            std::cerr << "\td^2/dydz\n";
            #endif
            const auto true_deriv = -sin(x)*sin(y)*cos(z);
            const auto calc_deriv = smart_deriv<0, 1, 1>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^2/dydz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }




        {   //d^3/dx^3
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dx^3\n";
            #endif
            const auto true_deriv = -cos(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<3, 0, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dx^3 failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dy^3
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dy^3\n";
            #endif
            const auto true_deriv = sin(x)*sin(y)*sin(z);
            const auto calc_deriv = smart_deriv<0, 3, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dy^3 failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dz^3
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dz^3\n";
            #endif
            const auto true_deriv = -sin(x)*cos(y)*cos(z);
            const auto calc_deriv = smart_deriv<0, 0, 3>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dz^3 failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

        {   //d^3/dx^2dy
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dx^2dy\n";
            #endif
            const auto true_deriv = sin(x)*sin(y)*sin(z);
            const auto calc_deriv = smart_deriv<2, 1, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dx^2dy failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dy^2dx
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dy^2dx\n";
            #endif
            const auto true_deriv = -cos(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<1, 2, 0>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dy^2dz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dx^2dz
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dx^2dz\n";
            #endif
            const auto true_deriv = -sin(x)*cos(y)*cos(z);
            const auto calc_deriv = smart_deriv<2, 0, 1>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dx^2dz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dz^2dx
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dz^2dx\n";
            #endif
            const auto true_deriv = -cos(x)*cos(y)*sin(z);
            const auto calc_deriv = smart_deriv<1, 0, 2>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dz^2dx failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dz^2dy
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dz^2dy\n";
            #endif
            const auto true_deriv = sin(x)*sin(y)*sin(z);
            const auto calc_deriv = smart_deriv<0, 1, 2>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dz^2dy failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dy^2dz
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dy^2dz\n";
            #endif
            const auto true_deriv = -sin(x)*cos(y)*cos(z);
            const auto calc_deriv = smart_deriv<0, 2, 1>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dy^2dz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
        {   //d^3/dxdydz
            #ifdef BVED_DLOG
            std::cerr << "\td^3/dxdydz\n";
            #endif
            const auto true_deriv = -cos(x)*sin(y)*cos(z);
            const auto calc_deriv = smart_deriv<1, 1, 1>(v, i).x();
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "d^3/dxdydz failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }





    }

    std::cout << "checking complete\n";
}

#ifdef BVED_DLOG
    #undef BVED_DLOG
#endif


#endif //CODE_BIG_VEC_DERIVS_HPP
