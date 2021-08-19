//
// Created by jacob on 8/8/21.
//

#ifndef CODE_DERIVATIVES_HPP
#define CODE_DERIVATIVES_HPP

#include "../MyMath/finite_difference.hpp"

namespace Examples {
    template <size_t x, size_t y, size_t z>
    struct test_d {
        double d[x][y][z];
        test_d(double d_[x][y][z]) {
            for (size_t i = 0; i < x; i++)
                for (size_t j = 0; j < y; j++)
                    for (size_t k = 0; k < z; k++) {
                        d[i][j][k] = d_[i][j][k];
                    }
        }

        double operator()(size_t i, size_t j, size_t k) const {return d[i][j][k];}
        double& operator()(size_t i, size_t j, size_t k) {return d[i][j][k];}
    };
}



void derivatives_ex() {
    constexpr unsigned dp = 20;
    constexpr double dx = 0.001;
    constexpr double dy = 0.0005;
    constexpr double dz = 0.002;
    constexpr double ex = 3;
    constexpr double ey = 4;
    constexpr double ez = 5;

    constexpr unsigned p = dp/2;

    double d[dp][dp][dp];

    for (size_t i = 0; i < dp; i++)
        for (size_t j = 0; j < dp; j++)
            for (size_t k = 0; k < dp; k++) {
                const auto x = static_cast<double>(i) - (dp)/2.0;
                const auto y = static_cast<double>(j) - (dp)/2.0;
                const auto z = static_cast<double>(k) - (dp)/2.0;

                d[i][j][k] = pow(sin(ex + x*dx),2) * pow(cos(ey + y*dy),2) * sin(ez + z*dz);
            }

    Examples::test_d<dp, dp, dp> data(d);


    //std::cout << data[p][p][p] << "\n";
    std::cout << "d/dx : \n";
    std::cout << "\tCentral \t" << central_difference_1st<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_1st<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_1st<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.114477\n";

    std::cout << "d/dy : \n";
    std::cout << "\tCentral \t" << central_difference_1st<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_1st<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_1st<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.0188936\n";

    std::cout << "d/dz : \n";
    std::cout << "\tCentral \t" << central_difference_1st<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_1st<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_1st<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.00241357\n";


    std::cout << "d^2/dx^2 : \n";
    std::cout << "\tCentral \t" << central_difference_2nd<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_2nd<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_2nd<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.786764\n";

    std::cout << "d^2/dy^2 : \n";
    std::cout << "\tCentral \t" << central_difference_2nd<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_2nd<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_2nd<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.00555718\n";

    std::cout << "d^2/dz^2 : \n";
    std::cout << "\tCentral \t" << central_difference_2nd<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_2nd<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_2nd<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.00815912\n";

    std::cout << "d^2/dxdy : \n";
    std::cout << "\tCentral \t" << central_difference_2nd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_2nd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_2nd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward(xy) \t" << central_forward_difference_2nd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward(yx) \t" << central_forward_difference_2nd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward(xy) \t" << central_backward_difference_2nd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward(yx) \t" << central_backward_difference_2nd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward(xy) \t" << forward_backward_difference_2nd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward(yx) \t" << forward_backward_difference_2nd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.265087\n";

    std::cout << "d^2/dxdz : \n";
    std::cout << "\tCentral \t" << central_difference_2nd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_2nd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_2nd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward(xz) \t" << central_forward_difference_2nd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward(zx) \t" << central_forward_difference_2nd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward(xz) \t" << central_backward_difference_2nd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward(zx) \t" << central_backward_difference_2nd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward(xz) \t" << forward_backward_difference_2nd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward(zx) \t" << forward_backward_difference_2nd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.0338637\n";

    std::cout << "d^2/dydz : \n";
    std::cout << "\tCentral \t" << central_difference_2nd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_2nd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_2nd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward(yz) \t" << central_forward_difference_2nd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward(zy) \t" << central_forward_difference_2nd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward(yz) \t" << central_backward_difference_2nd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward(zy) \t" << central_backward_difference_2nd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward(yz) \t" << forward_backward_difference_2nd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward(zy) \t" << forward_backward_difference_2nd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.00558898\n";


    std::cout << "d^3/dx^3 : \n";
    std::cout << "\tCentral \t" << central_difference_3rd<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd<0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.457907\n";

    std::cout << "d^3/dy^3 : \n";
    std::cout << "\tCentral \t" << central_difference_3rd<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd<1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.0755745\n";

    std::cout << "d^3/dz^3 : \n";
    std::cout << "\tCentral \t" << central_difference_3rd<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd<2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.00241357\n";



    std::cout << "d^3/dx^2dy : \n";
    std::cout << "\tCentral \t" << central_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward \t" << central_forward_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward \t" << central_backward_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-central \t" << forward_central_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward \t" << forward_backward_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-central \t" << backward_central_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-forward \t" << backward_forward_difference_3rd_mixed<0,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t1.82186\n";

    std::cout << "d^3/dx^2dz : \n";
    std::cout << "\tCentral \t" << central_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward \t" << central_forward_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward \t" << central_backward_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-central \t" << forward_central_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward \t" << forward_backward_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-central \t" << backward_central_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-forward \t" << backward_forward_difference_3rd_mixed<0,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.232735\n";

    std::cout << "d^3/dy^2dz : \n";
    std::cout << "\tCentral \t" << central_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward \t" << central_forward_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward \t" << central_backward_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-central \t" << forward_central_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward \t" << forward_backward_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-central \t" << backward_central_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-forward \t" << backward_forward_difference_3rd_mixed<1,2>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.00164389\n";

    std::cout << "d^3/dxdy^2 : \n";
    std::cout << "\tCentral \t" << central_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward \t" << central_forward_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward \t" << central_backward_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-central \t" << forward_central_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward \t" << forward_backward_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-central \t" << backward_central_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-forward \t" << backward_forward_difference_3rd_mixed<1,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t0.0779701\n";

    std::cout << "d^3/dxdz^2 : \n";
    std::cout << "\tCentral \t" << central_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward \t" << central_forward_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward \t" << central_backward_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-central \t" << forward_central_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward \t" << forward_backward_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-central \t" << backward_central_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-forward \t" << backward_forward_difference_3rd_mixed<2,0>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.114477\n";

    std::cout << "d^3/dydz^2 : \n";
    std::cout << "\tCentral \t" << central_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward \t" << forward_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward \t" << backward_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-forward \t" << central_forward_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tCentral-backward \t" << central_backward_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-central \t" << forward_central_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tForward-backward \t" << forward_backward_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-central \t" << backward_central_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tBackward-forward \t" << backward_forward_difference_3rd_mixed<2,1>(data, p, p, p, dx, dy, dz) << "\n";
    std::cout << "\tExpected\t-0.0188936\n";
}

#endif //CODE_DERIVATIVES_HPP
