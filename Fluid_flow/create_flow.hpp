//
// Created by jacob on 13/8/21.
//

#ifndef CODE_CREATE_FLOW_HPP
#define CODE_CREATE_FLOW_HPP


#include <fstream>
#include <chrono>

#include "make_mats.hpp"
#include "make_vecs.hpp"

#include "solver.hpp"
#include "../MyMath/boundary.hpp"

namespace flow {
    //width of the box
    constexpr double Wx = 3;
    constexpr double Wy = 4;
    constexpr double Wz = 5;

    //data points -1 along side of box
    constexpr unsigned N = 128;
    constexpr unsigned M = 128;
    constexpr unsigned P = 128;

    //size of grid
    constexpr double dx = Wx / static_cast<double>(N+1);
    constexpr double dy = Wy / static_cast<double>(M+1);
    constexpr double dz = Wz / static_cast<double>(P+1);

    //other consts
    constexpr double Re = 150; //http://www.airfoiltools.com/calculator/reynoldsnumber?MReNumForm%5Bvel%5D=10&MReNumForm%5Bchord%5D=0.2&MReNumForm%5Bkvisc%5D=1.3324E-5&yt0=Calculate
    constexpr double dt = 0.001;
    constexpr double max_t = 1;
}


template <unsigned N, unsigned M, unsigned P>
void set_BC(big_vec<N,M,P, vec3> &v, const double t) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            v.add_elm(i,j,0, 0,0,0);
            v.add_elm(i,j,P, 0,0,0);
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            v.add_elm(0,j,k, 0,0,0);
            v.add_elm(N,j,k, 0,0,0);
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            if (t < 0.2) {
                v.add_elm(i, 0, k, 0.1 * sin(0.1 * t), 0, 0);    //making the floor move
            } else {
                v.add_elm(i, 0, k, 0.1 * sin(0.1 * 0.2), 0, 0);
            }
            v.add_elm(i,M,k, 0,0,0);
        }
    }

}

template <unsigned N, unsigned M, unsigned P>
void enforce_PBC(big_vec<N,M,P, double> &p, const boundary_normals<N,M,P> &norms) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
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
                            p(i,j,k) = ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i+2,j,k)/3 + 4*p(i+1,j,k)/3;
                        } else if (!p.has_right(i,j,k)) {   //backward difference
                            p(i,j,k) =-ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i-2,j,k)/3 + 4*p(i-1,j,k)/3;
                        } else {    //central difference
                            p(i+1,j,k) = -ny/nx * smart_deriv<0,1,0>(p, i,j,k)*dx - nz/nx * smart_deriv<0,0,1>(p, i,j,k)*dx + p(i-1,j,k);
                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p.has_down(i, j, k)) { //forward difference
                            p(i,j,k) = nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j+2,k)/3 + 4*p(i,j+1,k)/3;
                        } else if (!p.has_up(i, j, k)) {  //backward difference
                            p(i,j,k) =-nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j-2,k)/3 + 4*p(i,j-1,k)/3;
                        } else {  //central difference
                            p(i,j+1,k) = -nx/ny * smart_deriv<1,0,0>(p, i,j,k)*dy - nz/ny * smart_deriv<0,0,1>(p, i,j,k)*dy + p(i,j-1,k);
                        }
                    } else {    //z direction
#ifndef NDEBUG
                        if (big_dir > 2) {
                            std::cerr << "the biggest direction cannot be larger than 2\n";
                        }
#endif
                        if (!p.has_front(i,j,k)) {  //forward difference
                            p(i,j,k) = ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k+2)/3 + 4*p(i,j,k+1)/3;
                        } else if (!p.has_back(i,j,k)) { //backward difference
                            p(i,j,k) =-ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k-2)/3 + 4*p(i,j,k-1)/3;
                        } else {
                            p(i,j,k+1) = -ny/nz * smart_deriv<0,1,0>(p, i,j,k)*dz - nx/nz * smart_deriv<1,0,0>(p, i,j,k)*dz + p(i,j,k-1);
                        }
                    }

                }

            }
        }
    }
}

template <unsigned N, unsigned M, unsigned P>
void v_IC(big_vec<N,M,P, vec3> &v) {
     //just leave it as 0
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                v.add_elm(i, j, k, 0, 0, 0);
            }
        }

    }
}

//num is the number of boundary points
template <unsigned N, unsigned M, unsigned P>
void create_boundary_points(boundary_points<N,M,P> &bound, unsigned &num) {
    num = 0;
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(i,0,k).has_down = false;
            bound(i,M,k).has_up = false;
            ++num;
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            bound(i,j,0).has_front = false;
            bound(i,j,P).has_back = false;
            ++num;
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(0,j,k).has_left = false;
            bound(N,j,k).has_right = false;
            ++num;
        }
    }
}


template <unsigned N, unsigned M, unsigned P>
void create_boundary_normals(boundary_normals<N,M,P> &bound) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            bound.add_point(i,0,k, vec3(0,1,0));
            bound.add_point(i,M,k, vec3(0,-1,0));
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            bound.add_point(i,j,0, vec3(0,0,1));
            bound.add_point(i,j,P, vec3(0,0,-1));
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            bound.add_point(0,j,k,vec3(1,0,0));
            bound.add_point(N,j,k,vec3(-1,0,0));
        }
    }
}


template <unsigned N, unsigned M, unsigned P>
void DEBUG_check_normal_for_all_boundary_points(const boundary_points<N,M,P> &point, const boundary_normals<N,M,P> &norm) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <=M; j++) {
            for (unsigned k = 0; k<=P; k++) {
                if (point.is_boundary(i,j,k)) {
                    if (!norm.contains(i,j,k)) {
                        std::cerr << "Boundary at i=" << i << " j=" << j << " k=" << k << "does not have a normal\n";
                    }
                }
            }
        }
    }
}


template<unsigned N, unsigned M, unsigned P, typename T>
void write_vec(const big_vec<N,M,P,T>& v, const char* file_loc) {
    std::ofstream output(file_loc);
    if (output.is_open()) {
        for (unsigned i = 0; i <= N; i++) {
            for (unsigned j = 0; j <= M; j++) {
                output << i * v.dx << " " << j * v.dy << " " << v(i, j, P / 2) << "\n";
            }
        }
    } else {
        std::cerr << "failed to open file\n";
    }

    output.close();
}

void solve_flow() {
    using namespace flow;

    std::cout << "Creating boundary";
    auto start = std::chrono::high_resolution_clock::now();
    unsigned no_boundary_points = 0;
    boundary_points<N,M,P> bound;
    create_boundary_points(bound, no_boundary_points);

    boundary_normals<N,M,P> norms(no_boundary_points);
    create_boundary_normals(norms);

#ifndef NDEBUG
    DEBUG_check_normal_for_all_boundary_points(bound, norms);
#endif


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Default constructing vectors";
    start = std::chrono::high_resolution_clock::now();
    big_vec<N,M,P,vec3> v_n(dx, dy, dz, &bound);
    big_vec<N,M,P,vec3> v_n1(dx, dy, dz, &bound);
    big_vec<N,M,P,double> p(dx, dy, dz, &bound);
    big_vec<N,M,P,vec3> bc(dx, dy, dz, &bound);
    big_vec<N,M,P,double> p_bc(dx, dy, dz, &bound);
    big_vec<N,M,P,double> p_c(dx, dy, dz, &bound);    //pressure correction vector


    big_vec<N,M,P,vec3> b(dx, dy, dz, &bound);
    big_vec<N,M,P,double> s(dx, dy, dz, &bound);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";



    std::cout << "Default constructing matrices";
    start = std::chrono::high_resolution_clock::now();
    big_matrix<N,M,P> Q(16);
    big_matrix<N,M,P> A(7);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Default constructing solvers";
    start = std::chrono::high_resolution_clock::now();

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";




    //initially creating matrices
    std::cout << "Creating A";
    start = std::chrono::high_resolution_clock::now();
    make_A(A, v_n, dt, Re);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Creating Q";
    start = std::chrono::high_resolution_clock::now();

    make_Q(Q, p, norms);

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Setting IC";
    start = std::chrono::high_resolution_clock::now();
    //first set IC
    v_IC(v_n);
    write_vec(v_n, "../DEBUG/velocity_data/0000.txt");
    write_vec(p, "../DEBUG/pressure_data/0000.txt");
    v_n1 = v_n;
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";


    std::cout << "Creating s";
    start = std::chrono::high_resolution_clock::now();
    //then create the s matrix
    enforce_PBC(p_bc, norms);
    make_s_first(s, Re, dt, v_n, p, p_bc);

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Solving for p";
    start = std::chrono::high_resolution_clock::now();
    //solve for p for the first timestep
    solve(Q, s, p_c);
    enforce_PBC(p_c, norms);
    p += p_c;

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    //setting BC vector
    set_BC(bc, 0);
    std::cout << "Making b";
    start = std::chrono::high_resolution_clock::now();
    //then make b
    make_b_first(b, Re, dt, v_n, p, bc);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Solving for v";
    start = std::chrono::high_resolution_clock::now();
    //and solve for v at the first timestep

    solve(A, b.xv, v_n.xv);
    solve(A, b.yv, v_n.yv);
    solve(A, b.zv, v_n.zv);

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";


    std::cout << "Enforcing BC";
    start = std::chrono::high_resolution_clock::now();
    //enforcing BC
    set_BC(v_n, 0);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";


    write_vec(v_n, "../DEBUG/velocity_data/0001.txt");
    write_vec(p, "../DEBUG/pressure_data/0001.txt");


    unsigned counter = 1;

    for (double t = dt; t < flow::max_t; t+=dt) {
        const auto start_loop = std::chrono::high_resolution_clock::now();
        const std::time_t start_time = std::chrono::system_clock::to_time_t(start_loop);
        std::cout << "t: " << t-dt << " / " << max_t << " at " << std::ctime(&start_time) << std::flush;

        std::cout << "\tMaking s" << std::flush;
        start = std::chrono::high_resolution_clock::now();
        //first make the s matrix
        enforce_PBC(p_bc, norms);
        make_s(s, Re, dt, v_n, v_n1, p, p_bc);

        //make_s_first(s, Re, dt, v_n);
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n" << std::flush;

        std::cout << "\tSolving for p" << std::flush;
        start = std::chrono::high_resolution_clock::now();
        //solve for p for the next timestep

        solve(Q, s, p_c);
        enforce_PBC(p_c, norms);
        p += p_c;

        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n" << std::flush;

        //setting BC vector
        set_BC(bc, t);
        std::cout << "\tMaking b";
        start = std::chrono::high_resolution_clock::now();
        //then make b
        make_b(b, Re, dt, v_n, v_n1, p, bc);
        //make_b_first(b, Re, dt, v_n, p);
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        std::cout << "\tShuffling variables";
        start = std::chrono::high_resolution_clock::now();
        //shuffling of variables
        v_n1 = v_n;
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        //solve for v for the next time step
        std::cout << "\tSolving for vx" << std::flush;
        start = std::chrono::high_resolution_clock::now();
        solve(A, b.xv, v_n.xv);

        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";


        std::cout << "\tSolving for vy" << std::flush;
        start = std::chrono::high_resolution_clock::now();
        solve(A, b.yv, v_n.yv);

        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        std::cout << "\tSolving for vz" << std::flush;
        start = std::chrono::high_resolution_clock::now();
        solve(A, b.zv, v_n.zv);

        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";


        std::cout << "\tEnforcing BC";
        start = std::chrono::high_resolution_clock::now();
        //enforcing BC
        set_BC(v_n, t);
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        ++counter;
        std::string file_name;
        if (counter < 10) {
            file_name = "000" + std::to_string(counter);
        } else if (counter < 100) {
            file_name = "00" + std::to_string(counter);
        } else if (counter < 1000) {
            file_name = "0" + std::to_string(counter);
        } else {
            file_name = std::to_string(counter);
        }

        std::cout << "\tWriting result " << file_name;
        start = std::chrono::high_resolution_clock::now();
        write_vec(v_n, ("../DEBUG/velocity_data/" + file_name + ".txt").c_str() );
        write_vec(p, ("../DEBUG/pressure_data/" + file_name + ".txt").c_str() );
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        const auto end_loop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> dur_loop = end_loop - start_loop;
        std::cout << "Timestep took : " << dur_loop.count() << "s\n";
    }


}



#endif //CODE_CREATE_FLOW_HPP
