//
// Created by jacob on 13/8/21.
//

#ifndef CODE_CREATE_FLOW_HPP
#define CODE_CREATE_FLOW_HPP


#include <fstream>
#include <chrono>
#include <string>
#include <ctime>
#include <ratio>

#include "make_mats.hpp"
#include "make_vecs.hpp"

#include "solver.hpp"
#include "../MyMath/boundary.hpp"
#include "flow_env.hpp"
#include "timing.hpp"


struct output_settings {
    std::string_view vel_file_loc = "../DEBUG/velocity_data/";
    std::string_view pres_file_loc = "../DEBUG/pressure_data/";
    std::string_view time_file_name = "../DEBUG/times.txt";

      //only used if write_all_times is false
    std::string_view final_vel_name = "../DEBUG/vel_final.txt";
    std::string_view final_pres_name = "../DEBUG/vel_final.txt";
};

//for choice of Reynolds number see //http://www.airfoiltools.com/calculator/reynoldsnumber?MReNumForm%5Bvel%5D=10&MReNumForm%5Bchord%5D=0.2&MReNumForm%5Bkvisc%5D=1.3324E-5&yt0=Calculate
//Wx,Wy,Wz represent the width of the box
template<unsigned no_timesteps, unsigned N, unsigned M, unsigned P, bool write_all_times = true>
void solve_flow(const output_settings &os, const double max_t = 1, const double Re = 150, const double Wx = 3, const double Wy = 4, const double Wz = 5) {
    //size of grid
    double dx = Wx / static_cast<double>(N+1);
    double dy = Wy / static_cast<double>(M+1);
    double dz = Wz / static_cast<double>(P+1);

    double dt = max_t / (double)no_timesteps;

    flow_timer timer(os.time_file_name.data() );


    unsigned no_boundary_points = 0;
    boundary_points<N,M,P> bound;
    create_boundary_points(bound, no_boundary_points);

    boundary_normals<N,M,P> norms(no_boundary_points);
    create_boundary_normals(norms);

#ifndef NDEBUG
    DEBUG_check_normal_for_all_boundary_points(bound, norms);
#endif


    big_vec<N,M,P,vec3> v_n(dx, dy, dz, &bound);    //velocity at hte current time-step
    big_vec<N,M,P,vec3> v_n1(dx, dy, dz, &bound);   //velocity at the previous time-step
    big_vec<N,M,P,double> p(dx, dy, dz, &bound);    //pressure vector
    big_vec<N,M,P,vec3> bc(dx, dy, dz, &bound);     //boundary conditions (used to set b)
    big_vec<N,M,P,double> p_bc(dx, dy, dz, &bound); //pressure boundary conditions (used to set s)
    big_vec<N,M,P,double> p_c(dx, dy, dz, &bound);    //pressure correction vector


    big_vec<N,M,P,vec3> b(dx, dy, dz, &bound);
    big_vec<N,M,P,double> s(dx, dy, dz, &bound);


    big_matrix<N,M,P> Q(16);
    big_matrix<N,M,P> A(7);


    //initially creating matrices
    make_A(A, v_n, dt, Re);

    //making Q
    make_Q(Q, p, norms);


    //first set IC
    v_IC(v_n);
    if constexpr (write_all_times) {
        write_vec(v_n, (std::string(os.vel_file_loc) + "0000.txt").data());
        write_vec(p, (std::string(os.pres_file_loc) + "0000.txt").data());
    }
    v_n1 = v_n;

    //then create the s matrix
    enforce_PBC(p_bc, norms);
    make_s_first(s, Re, dt, v_n, p, p_bc);

    //solve for p for the first timestep
    solve(Q, s, p_c);
    enforce_PBC(p_c, norms);
    p += p_c;

    //setting BC vector
    set_BC(bc, 0);
    //then make b
    make_b_first(b, Re, dt, v_n, p, bc);

    //and solve for v at the first timestep

    solve(A, b.xv, v_n.xv);
    solve(A, b.yv, v_n.yv);
    solve(A, b.zv, v_n.zv);

    //enforcing BC
    set_BC(v_n, 0);

    if constexpr (write_all_times) {
        write_vec(v_n, (std::string(os.vel_file_loc) + "0001.txt").data());
        write_vec(p, (std::string(os.pres_file_loc) + "0001.txt").data());
    }


    unsigned counter = 1;

    for (double t = dt; t < max_t; t+=dt) {
        const auto start_loop = std::chrono::high_resolution_clock::now();
        const std::time_t start_time = std::chrono::system_clock::to_time_t(start_loop);
        char time_human[9]; //9 characters for HH:MM:SS (8char) plus termination \0
        if (std::strftime(time_human, sizeof(time_human), "%T", std::localtime(&start_time))) {
            std::cout << "t: " << t << " / " << max_t << " at " << time_human << std::flush;
        } else {
            std::cerr << "Timing error\n";
        }


        timer.set_start(std::chrono::high_resolution_clock::now());
        //first make the s matrix
        enforce_PBC(p_bc, norms);
        make_s(s, Re, dt, v_n, v_n1, p, p_bc);

        timer.set_end(std::chrono::high_resolution_clock::now());
        timer.save_s_create_time();


        timer.set_start(std::chrono::high_resolution_clock::now());
        //solve for p for the next timestep
        solve(Q, s, p_c);
        enforce_PBC(p_c, norms);
        p += p_c;
        timer.set_end(std::chrono::high_resolution_clock::now());
        timer.save_p_solve_time();

        //setting BC vector
        set_BC(bc, t);

        timer.set_start(std::chrono::high_resolution_clock::now());
        //then make b
        make_b(b, Re, dt, v_n, v_n1, p, bc);
        timer.set_end(std::chrono::high_resolution_clock::now());
        timer.save_b_create_time();


        //shuffling of variables
        v_n1 = v_n;


        //solve for v for the next time step
        timer.set_start(std::chrono::high_resolution_clock::now());
        solve(A, b.xv, v_n.xv);
        timer.set_end(std::chrono::high_resolution_clock::now());
        timer.save_vx_solve_time();


        timer.set_start(std::chrono::high_resolution_clock::now());
        solve(A, b.yv, v_n.yv);
        timer.set_end(std::chrono::high_resolution_clock::now());
        timer.save_vy_solve_time();

        timer.set_start(std::chrono::high_resolution_clock::now());
        solve(A, b.zv, v_n.zv);
        timer.set_end(std::chrono::high_resolution_clock::now());
        timer.save_vz_solve_time();


        //enforcing BC
        set_BC(v_n, t);


        ++counter;
        std::string file_name;
        if constexpr (write_all_times) {
            if (counter < 10) {
                file_name = "000" + std::to_string(counter);
            } else if (counter < 100) {
                file_name = "00" + std::to_string(counter);
            } else if (counter < 1000) {
                file_name = "0" + std::to_string(counter);
            } else {
                file_name = std::to_string(counter);
            }


            write_vec(v_n, (std::string(os.vel_file_loc) + file_name + ".txt").data());
            write_vec(p, (std::string(os.pres_file_loc) + file_name + ".txt").data());
        }

        const auto end_loop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> dur_loop = end_loop - start_loop;
        std::cout << "\ttimestep took : " << dur_loop.count() << "s";

        if constexpr (write_all_times) {
            std::cout << "\tfile written : " << file_name << "\n";
        } else {
            std::cout << "\n";
        }

        timer.write_times(t);
    }

    if constexpr (!write_all_times) {
        write_vec(v_n, os.final_vel_name.data());
        write_vec(p, os.final_pres_name.data());
    }


}



#endif //CODE_CREATE_FLOW_HPP
