//
// Created by jacob on 13/8/21.
//

#ifndef CODE_CREATE_FLOW_HPP
#define CODE_CREATE_FLOW_HPP


#include <fstream>
#include <chrono>
#include <string>

#include "make_mats.hpp"
#include "make_vecs.hpp"

#include "solver.hpp"
#include "../MyMath/boundary.hpp"
#include "flow_env.hpp"

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

    constexpr std::string_view vel_file_loc = "../DEBUG/velocity_data/";
    constexpr std::string_view pres_file_loc = "../DEBUG/pressure_data/";
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
    write_vec(v_n, (std::string(vel_file_loc) + "0000.txt").data());
    write_vec(p, (std::string(pres_file_loc) + "0000.txt").data());
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


    write_vec(v_n, (std::string(vel_file_loc) + "0001.txt").data());
    write_vec(p, (std::string(pres_file_loc) + "0001.txt").data());


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
        write_vec(v_n, (std::string(vel_file_loc) + file_name + ".txt").c_str() );
        write_vec(p, (std::string(pres_file_loc) + file_name + ".txt").c_str() );
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        const auto end_loop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> dur_loop = end_loop - start_loop;
        std::cout << "Timestep took : " << dur_loop.count() << "s\n";
    }


}



#endif //CODE_CREATE_FLOW_HPP
