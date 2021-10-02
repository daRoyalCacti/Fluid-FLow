//
// Created by jacob on 3/10/21.
//

#ifndef CODE_INTERP_HPP
#define CODE_INTERP_HPP

#include "../MyMath/calc.hpp"
#include "../MyMath/boundary.hpp"

void interp_ex() {



    constexpr unsigned nx = 20;
    constexpr unsigned ny = 26;
    constexpr unsigned nz = 22;

    constexpr double dx = 0.001;
    constexpr double dy = 0.0004;
    constexpr double dz = 0.0007;

    constexpr double x_off = dx/2;
    constexpr double y_off = dy/3;
    constexpr double z_off = 2*dz/3;

    constexpr double ex = 3;
    constexpr double ey = 4;
    constexpr double ez = 5;


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

    std::cout << "testing to see if interpolation correct\n";

    boundary_conditions BC(&m, dx, dy, dz, nx*dx, ny*dy, nz*dz, ex, ey, ez);
    BC.global_grid.DEBUG_write_boundary_points();
    BC.DEBUG_write_normal_vectors();
    //BC.global_grid.DEBUG_write_boundary_points_at_x(3.007);
    //BC.DEBUG_write_normal_vectors_at_x(3.007);
    big_vec_d v(BC);

    for (unsigned i = 0; i < v.size(); i++) {
        const auto x = BC.global_grid.x[i];
        const auto y = BC.global_grid.y[i];
        const auto z = BC.global_grid.z[i];
        v(i) = cos(x)*sin(y)*sin(z);
    }

    v.move(x_off, y_off, z_off);

    constexpr double tol = 0.1;     //make smaller
    for (unsigned i = 0; i < v.size(); i++) {
        const auto x = BC.global_grid.x[i] + x_off;
        const auto y = BC.global_grid.y[i] + y_off;
        const auto z = BC.global_grid.z[i] + z_off;

        if (!std::isfinite(v(i)) || std::abs(v(i) - cos(x)*sin(y)*sin(z)) > tol ) {
            std::cerr << "interpolation failed at ind = " << i << "\n";
            std::cerr << "\ttrue value = " << cos(x)*sin(y)*sin(z) << "\tinterpolated value = " << v(i) << "\n";
            std::cerr << "\tind " << i << " can move: ";
            if (v.can_move(i,-1,0,0)) std::cerr << " left, ";
            if (v.can_move(i,1,0,0)) std::cerr << " right, ";
            if (v.can_move(i,0,-1,0)) std::cerr << " down, ";
            if (v.can_move(i,0,1,0)) std::cerr << " up, ";
            if (v.can_move(i,0,0,-1)) std::cerr << " forward, ";
            if (v.can_move(i,0,0,1)) std::cerr << " backward, ";
            std::cerr << "\n";
        }
    }

    std::cout << "testing complete\n";
}

#endif //CODE_INTERP_HPP
