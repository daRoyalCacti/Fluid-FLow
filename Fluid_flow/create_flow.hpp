//
// Created by jacob on 13/8/21.
//

#ifndef CODE_CREATE_FLOW_HPP
#define CODE_CREATE_FLOW_HPP

#define VC
//#define NO_MESH_UPDATE
#define FLUID_MOVES_MESH
#define SAME_PLOTTING_INDS
#define DLOG    //detailed logging
//#define OFFSET_MESH_UPDATE  //if the mesh should only be updated every 2nd timestep
//#define REPEATS
#ifdef REPEATS
constexpr unsigned no_repeats = 3;
#endif

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
#include "boundary_conditions.hpp"
#include "update_vecs.hpp"

#include "create_grids.hpp"
#include "../Rigid_body/body.hpp"
#include "../Rigid_body/triangle_mesh.hpp"
#include "../MyMath/boundary.hpp"
#include "../MyMath/big_vec.hpp"
#include "update_mesh.hpp"


constexpr double accuracy_percent_v = 0.0001;
constexpr double accuracy_percent_p = 0.01;

struct output_settings {
    std::string_view vel_file_loc = "../DEBUG/velocity_data/";
    std::string_view pres_file_loc = "../DEBUG/pressure_data/";
    std::string_view time_file_name = "../DEBUG/times.txt";
    std::string_view body_file_loc = "../DEBUG/rigid_body_data/";

      //only used if write_all_times is false
    std::string_view final_vel_name = "../DEBUG/vel_final.txt";
    std::string_view final_pres_name = "../DEBUG/vel_final.txt";
};


//for choice of Reynolds number see //http://www.airfoiltools.com/calculator/reynoldsnumber?MReNumForm%5Bvel%5D=10&MReNumForm%5Bchord%5D=0.2&MReNumForm%5Bkvisc%5D=1.3324E-5&yt0=Calculate
//Wx,Wy,Wz represent the width of the box
template<unsigned no_timesteps, unsigned N, unsigned M, unsigned P, bool write_all_times = true>
void solve_flow(body *rb, const output_settings &os, const double max_t = 1, const double Re = 150, const double Wx = 3, const double Wy = 4, const double Wz = 5) {
    //size of grid
    double dx = Wx / static_cast<double>(N+1);
    double dy = Wy / static_cast<double>(M+1);
    double dz = Wz / static_cast<double>(P+1);

    double dt = max_t / (double)no_timesteps;

    //creating timer
    flow_timer timer(os.time_file_name.data() );



#ifdef DLOG
    std::cout << "setting boundary conditions\n";
#endif
    boundary_conditions BC(&rb->model, dx, dy, dz, Wx, Wy, Wz);
    BC.global_grid.DEBUG_write_boundary_points(false);
    BC.DEBUG_write_normal_vectors();



 #ifdef DLOG
 std::cout << "constructing matrices and vectors\n";
 #endif
     big_vec_v v_n(BC);    //velocity at the current time-step
     big_vec_v v_n1(BC);   //velocity at the previous time-step
#ifdef VC
    big_vec_v vc(BC);
#endif

     big_vec_d p(BC);    //pressure vector
     big_vec_d p_c(BC);    //pressure correction vector



     big_vec_v b(BC);
     big_vec_d s(BC);


     big_matrix Q(BC, 16);
     big_matrix A(BC, 7);
     //return;
 #ifdef DLOG
     std::cout << "creating A\n";
 #endif
     //initially creating matrices
     make_A(A, v_n, dt, Re, BC);

 #ifdef DLOG
     std::cout << "creating Q\n";
 #endif
     //making Q
     make_Q(Q, p, BC);

 #ifdef DLOG
     std::cout << "setting velocity IC\n";
 #endif
     //first set IC
     v_IC(v_n);

#ifdef SAME_PLOTTING_INDS
     const auto inds = v_n.g->get_middle_inds();
#endif


#ifdef DLOG
     std::cout << "writing velocity IC\n";
#endif
     if constexpr (write_all_times) {
#ifdef DLOG
         std::cout << "\tGetting indices\n";
#endif


#ifndef SAME_PLOTTING_INDS
         const auto inds = v_n.g->get_middle_inds();
#endif


#ifdef DLOG
         std::cout << "\tWriting v\n";
 #endif
         write_vec(v_n,inds, (std::string(os.vel_file_loc) + "0000.txt").data());
 #ifdef DLOG
         std::cout << "\tWriting p\n";
 #endif
         write_vec(p, inds,(std::string(os.pres_file_loc) + "0000.txt").data());
 #ifdef DLOG
         std::cout << "\tWriting r\n";
 #endif
         rb->write_pos((std::string(os.body_file_loc) + "0000.txt").data());
     }
#ifdef DLOG
    std::cout << "copying vector\n";
#endif
     v_n1 = v_n;

#ifndef     NO_MESH_UPDATE
#ifdef STORE_SOLVERS
    interp_solvers interps(BC.global_grid);
    const grid init_grid = BC.global_grid;  //creating a copy of the grid that doesn't update for the interpolation
    const vec3 init_com = rb->model.pos_cm;
#endif

 #ifdef DLOG
     std::cout << "setting the mesh boundary conditions\n";
 #endif
#ifdef STORE_SOLVERS
     update_mesh(BC, rb, v_n, v_n1, p, dt, 0, interps, init_grid, init_com);
#else
     update_mesh(BC, rb, v_n, v_n1, p, dt, 0);
#endif
#endif

 #ifdef DLOG
     std::cout << "making s\n";
 #endif
     make_s_first(s, Re, dt, v_n, p);

 #ifdef DLOG
     std::cout << "Solving pressure\n";
 #endif
     //solve for p for the first timestep
     solve_pressure(BC, Q, s, p_c, p, accuracy_percent_p);


     //setting BC vector
 #ifdef DLOG
     std::cout << "making b\n";
 #endif
     //then make b
     make_b_first(b, Re, dt, v_n, p, BC);

 #ifdef DLOG
     std::cout << "solving for v\n";
 #endif
     //and solve for v at the first timestep
     solve_velocity(BC, A, b, vc, v_n, accuracy_percent_v);

 #ifdef DLOG
     std::cout << "writing first timestep\n";
 #endif
     if constexpr (write_all_times) {
 #ifdef DLOG
         std::cout << "\tGetting indices\n";
 #endif

#ifndef SAME_PLOTTING_INDS
         const auto inds = v_n.g->get_middle_inds();
#endif


#ifdef DLOG
         std::cout << "\tWriting v\n";
 #endif
         write_vec(v_n,inds, (std::string(os.vel_file_loc) + "0001.txt").data());
 #ifdef DLOG
         std::cout << "\tWriting p\n";
 #endif
         write_vec(p, inds,(std::string(os.pres_file_loc) + "0001.txt").data());
 #ifdef DLOG
         std::cout << "\tWriting r\n";
 #endif
         rb->write_pos((std::string(os.body_file_loc) + "0001.txt").data());
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

#ifndef NO_MESH_UPDATE
         //updating the mesh
         timer.set_start(std::chrono::high_resolution_clock::now());
#ifdef OFFSET_MESH_UPDATE
         if (counter % 2 == 0) {
#endif

             //only updating mesh every 2nd frame
             // - solve fluid at old position
             // - update mesh
             // - solve fluid at new position
             // - repeat

#ifdef STORE_SOLVERS
            update_mesh(BC, rb, v_n, v_n1, p, dt, t, interps, init_grid, init_com);
#else
             update_mesh(BC, rb, v_n, v_n1, p, dt, t);
#endif


#ifdef OFFSET_MESH_UPDATE
         }
#endif
             timer.set_end(std::chrono::high_resolution_clock::now());
             timer.save_mesh_update_time();
#endif


#ifdef REPEATS
         for (unsigned i = 0; i < no_repeats; i++) {
#endif
             //for (unsigned j = 0; j < 3; j++) {
                 timer.set_start(std::chrono::high_resolution_clock::now());
                 //first make the s matrix
                 make_s(s, Re, dt, v_n, v_n1, p);

                 timer.set_end(std::chrono::high_resolution_clock::now());
                 timer.save_s_create_time();


                 timer.set_start(std::chrono::high_resolution_clock::now());
                 //solve for p for the next timestep

                 solve_pressure(BC, Q, s, p_c, p, accuracy_percent_p);
             //}

             timer.set_end(std::chrono::high_resolution_clock::now());
             timer.save_p_solve_time();


             timer.set_start(std::chrono::high_resolution_clock::now());
             //then make b
             make_b(b, Re, dt, v_n, v_n1, p, BC);
             timer.set_end(std::chrono::high_resolution_clock::now());
             timer.save_b_create_time();


             //shuffling of variables
             v_n1 = v_n;

             //solve for v for the next time step
             timer.set_start(std::chrono::high_resolution_clock::now());
             solve_velocity(BC, A, b, vc, v_n, accuracy_percent_v);
             timer.set_end(std::chrono::high_resolution_clock::now());
             timer.save_v_solve_time();


 #ifdef REPEATS
         }
 #endif


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
#ifndef SAME_PLOTTING_INDS
             const auto inds = v_n.g->get_middle_inds();
#endif
             write_vec(v_n,inds, (std::string(os.vel_file_loc) + file_name + ".txt").data());
             write_vec(p,inds, (std::string(os.pres_file_loc) + file_name + ".txt").data());
             rb->write_pos((std::string(os.body_file_loc) + file_name + ".txt").data());
         }

         const auto end_loop = std::chrono::high_resolution_clock::now();
         const std::chrono::duration<double> dur_loop = end_loop - start_loop;
         std::cout << "\ttimestep took : " << dur_loop.count() << "s";

         timer.save_total_time(start_loop, end_loop);

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
