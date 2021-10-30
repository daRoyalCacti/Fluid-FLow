//
// Created by jacob on 12/8/21.
//

#ifndef CODE_CALC_DERIVS_HPP
#define CODE_CALC_DERIVS_HPP

#include "../MyMath/calc.hpp"
#include "../MyMath/boundary.hpp"

#define BCD_DLOG   //bvec example derivating detailed logging

void calc_derivs_ex() {

    constexpr unsigned nx = 20;
    constexpr unsigned ny = 26;
    constexpr unsigned nz = 22;

    constexpr double dx = 0.001;
    constexpr double dy = 0.0004;
    constexpr double dz = 0.0007;

    constexpr double ex = 3;
    constexpr double ey = 4;
    constexpr double ez = 5;

    constexpr double tol = 0.1;    //tolerance for if the derivative is correct


    constexpr double mass1 = 0.01;

    constexpr vec3 vel_cm = vec3(0.75, 0, 0);
    constexpr vec3 w{};

    constexpr double wx = nx*dx;
    constexpr double wy = ny*dy;
    constexpr double wz = nz*dz;

    constexpr double z_mid = ez + wz/2;
    constexpr double y_mid = ey + wy/2;
    constexpr double x_mid = ex + wx/2;

    std::vector<vec3> pos;
    std::vector<unsigned> inds;
    std::vector<double> mass;
    std::vector<vec3> vels;
    std::vector<vec3> norms;

    const unsigned box_reduction = 10;

    //left/right rectangles
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                pos.emplace_back( (2*i-1)*wx/box_reduction + x_mid, (2*j-1)*wy/box_reduction + y_mid, (2*k-1)*wz/box_reduction + z_mid );
                mass.push_back(mass1);
                norms.emplace_back((2*i-1), 0, 0);
            }
        }
    }

    //top/bottom rectangles
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                pos.emplace_back( (2*i-1)*wx/box_reduction + x_mid, (2*j-1)*wy/box_reduction + y_mid, (2*k-1)*wz/box_reduction + z_mid );
                mass.push_back(mass1);
                norms.emplace_back(0, (2*j-1), 0);
            }
        }
    }

    //front/back rectangles
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
                pos.emplace_back( (2*i-1)*wx/box_reduction + x_mid, (2*j-1)*wy/box_reduction + y_mid, (2*k-1)*wz/box_reduction + z_mid );
                mass.push_back(mass1);
                norms.emplace_back(0, 0, (2*k-1));
            }
        }
    }

    for (int i = 0; i < 6; i++) {
        inds.push_back(4*i);    //lower left
        inds.push_back(4*i+1);  //lower right
        inds.push_back(4*i+2);  //upper left

        inds.push_back(4*i+1);  //lower right
        inds.push_back(4*i+2);  //upper left
        inds.push_back(4*i+3);  //upper right
    }

    //mesh m(pos, inds, mass, norms, vec3(0.75,0,0), vec3(0,0,0));
    mesh m(pos, inds, mass, norms, vec3(0.0,0,0), vec3(0,0,0));
    //body b(m);

    std::cout << "testing to see if derivatives correct\n";

    boundary_conditions BC(&m, dx, dy, dz, nx*dx, ny*dy, nz*dz, ex, ey, ez);
    BC.global_grid.DEBUG_write_boundary_points();
    BC.DEBUG_write_normal_vectors();
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
#ifdef BCD_DLOG
        std::cerr << i+1 << "/" << v.size() << "  \t" << x << " " << y << " " << z << "\n";
#endif

#ifdef BCD_DLOG
        std::cerr << "\tAdvection\n";
#endif
        {   //advection
            const auto true_deriv = vec3(cos(x) * pow(cos(z),2) *sin(x) *pow(sin(y),2) + cos(x) *pow(cos(y),2)* sin(x) *pow(sin(z),2) - cos(x) *sin(x)* pow(sin(y),2) *pow(sin(z),2),
                                         cos(y)* pow(cos(z),2) *pow(sin(x),2)* sin(y) + pow(cos(x),2)* cos(y) *sin(y) *pow(sin(z),2) - cos(y) *pow(sin(x),2) *sin(y) *pow(sin(z),2),
                                         pow(cos(y),2) *cos(z)* pow(sin(x),2) *sin(z) + pow(cos(x),2) *cos(z) *pow(sin(y),2) *sin(z) - cos(z) *pow(sin(x),2)* pow(sin(y),2) *sin(z));
            const auto calc_deriv = advection_old(v, i);
            if ((calc_deriv - true_deriv).length_squared() > tol) {
                std::cerr << "advection failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
#ifdef BCD_DLOG
        std::cerr << "\tlapacian\n";
#endif
        {   //lapacian
            const auto true_deriv = vec3(-3*cos(x)*sin(y)*sin(z), -3*cos(y)*sin(x)*sin(z), -3*cos(z)*sin(x)*sin(y));
            const auto calc_deriv = laplacian_old(v, i);
            if ((calc_deriv - true_deriv).length_squared() > tol) {
                std::cerr << "lapacian failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
#ifdef BCD_DLOG
        std::cerr << "\tDivergence\n";
#endif
        {   //Divergence
            const auto true_deriv =-3*sin(x)*sin(y)*sin(z);
            const auto calc_deriv = divergence_old(v, i);
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "Divergence failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
#ifdef BCD_DLOG
        std::cerr << "\tDivergence of advection\n";
#endif
        {   //Divergence of advection
            const auto true_deriv =2* pow(cos(y),2)* pow(cos(z),2) *pow(sin(x),2) + 2 *pow(cos(x),2) *pow(cos(z),2)* pow(sin(y),2) -
                    3* pow(cos(z),2) *pow(sin(x),2) *pow(sin(y),2) + 2 *pow(cos(x),2) *pow(cos(y),2) *pow(sin(z),2) -
                    3 *pow(cos(y),2)* pow(sin(x),2)* pow(sin(z),2) - 3 *pow(cos(x),2) *pow(sin(y),2) *pow(sin(z),2) +
                    3*pow(sin(x),2)* pow(sin(y),2) *pow(sin(z),2);
            const auto calc_deriv = divergence_advection_old(v, i);
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "Divergence of advection failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }
#ifdef BCD_DLOG
        std::cerr << "\tDivergence of laplacian\n";
#endif
        {   //Divergence of laplacian
            const auto true_deriv = 9 * sin(x) * sin(y) * sin(z);
            const auto calc_deriv = divergence_laplacian_old(v, i);
            if (abs(calc_deriv - true_deriv) > tol) {
                std::cerr << "Divergence of laplacian failed at ind = " << i << "\n";
                std::cerr << "true value : " << true_deriv << "\t calculated value : " << calc_deriv << "\n";
            }
        }

    }


    std::cout << "testing complete\n";
}

#ifdef BCD_DLOG
    #undef BCD_DLOG
#endif

#endif //CODE_CALC_DERIVS_HPP
