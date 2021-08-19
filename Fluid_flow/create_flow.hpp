//
// Created by jacob on 13/8/21.
//

#ifndef CODE_CREATE_FLOW_HPP
#define CODE_CREATE_FLOW_HPP

//#define VIENNACL_WITH_CUDA
//#define GPU_LIN_ALG
#define testing_gpu
//#define use_pc  //use pressure correction

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

    //size of gride
    constexpr double dx = Wx / static_cast<double>(N+1);
    constexpr double dy = Wy / static_cast<double>(M+1);
    constexpr double dz = Wz / static_cast<double>(P+1);

    //other consts
    constexpr double Re = 150; //http://www.airfoiltools.com/calculator/reynoldsnumber?MReNumForm%5Bvel%5D=10&MReNumForm%5Bchord%5D=0.2&MReNumForm%5Bkvisc%5D=1.3324E-5&yt0=Calculate
    constexpr double dt = 0.001;
    constexpr double max_t = 1;
}


template <unsigned N_, unsigned M_, unsigned P_>
void set_BC(big_vec<N_,M_,P_, vec3> &v, const double t) {
    for (unsigned i = 0; i <= N_; i++) {
        for (unsigned k = 0; k <= P_; k++) {
            v.add_elm(i,0,k, 0.1*sin(0.1*t), 0, 0);    //making the floor move
            v.add_elm(i,M_,k, 0,0,0);
        }
    }
    for (unsigned i = 0; i <= N_; i++) {
        for (unsigned j = 0; j <= M_; j++) {
            v.add_elm(i,j,0, 0,0,0);
            v.add_elm(i,j,P_, 0,0,0);
        }
    }
    for (unsigned j = 0; j <= M_; j++) {
        for (unsigned k = 0; k <= P_; k++) {
            v.add_elm(0,j,k, 0,0,0);
            v.add_elm(N_,j,k, 0,0,0);
        }
    }

}

template <unsigned N, unsigned M, unsigned P>
void v_IC(big_vec<N,M,P, vec3> &v) {
    //set_BC(v);  //very boring initial conditon --- currently just making it so the flow at the walls is set properly
    /*for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                v.add_elm(i,j,k, flow::wall_v * (1 - j / static_cast<double>(M) )   //linearly interpolating between moving wall and not moving wall
                );
            }
        }
    }*/

    //just leave it as 0
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                v.add_elm(i, j, k, 0, 0, 0);
            }
        }

    }
}

template <unsigned N, unsigned M, unsigned P>
void create_boundary_points(boundary_points<N,M,P> &bound) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(i,0,k).has_down = false;
            bound(i,M,k).has_up = false;
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            bound(i,j,0).has_front = false;
            bound(i,j,P).has_back = false;
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(0,j,k).has_left = false;
            bound(N,j,k).has_right = false;
        }
    }
}

//delete start
#include <string_view>

template <typename T>
constexpr auto type_name() noexcept {
    std::string_view name = "Error: unsupported compiler", prefix, suffix;
#ifdef __clang__
    name = __PRETTY_FUNCTION__;
  prefix = "auto type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
    name = __PRETTY_FUNCTION__;
    prefix = "constexpr auto type_name() [with T = ";
    suffix = "]";
#elif defined(_MSC_VER)
    name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void) noexcept";
#endif
    name.remove_prefix(prefix.size());
    name.remove_suffix(suffix.size());
    return name;
}
//delete end


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

    boundary_points<N,M,P> bound;
    create_boundary_points(bound);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Default constructing vectors";
    start = std::chrono::high_resolution_clock::now();
    big_vec<N,M,P,vec3> v_n(dx, dy, dz, &bound);
    big_vec<N,M,P,vec3> v_n1(dx, dy, dz, &bound);
    big_vec<N,M,P,double> p(dx, dy, dz, &bound);
    big_vec<N,M,P,vec3> bc(dx, dy, dz, &bound);
#ifdef use_pc
    big_vec<N,M,P,double> p_c(dx, dy, dz);    //pressure correction vector
#endif


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
#ifdef GPU_LIN_ALG
#ifdef testing_gpu

#else
    solver A_solver(A);
    solver Q_solver(Q);
#endif
#else
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::IncompleteLUT<double> > A_solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::IncompleteLUT<double>  > Q_solver;
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor> > A_solver;
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>  > Q_solver;
    //A_solver.setTolerance(0.001);
    //Q_solver.setTolerance(0.001);

    //A_solver.setMaxIterations(200);
    //Q_solver.setMaxIterations(200);
#endif
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
    make_Q(Q, p);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    //setting the solvers
#ifdef GPU_LIN_ALG

#else
    std::cout << "Setting the solvers";
    start = std::chrono::high_resolution_clock::now();
    A_solver.compute(A.m);
    Q_solver.compute(Q.m);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";
#endif


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
    make_s_first(s, Re, dt, v_n, p);
    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    std::cout << "\t" << dur.count() << "s\n";

    std::cout << "Solving for p";
    start = std::chrono::high_resolution_clock::now();
    //solve for p for the first timestep
    //p.v = Q_solver.solve(s.v);
#ifdef GPU_LIN_ALG
    #ifdef testing_gpu
        #ifdef use_pc
            solve(Q, s, p_c);
            p += p_c;
        #else
            solve(Q, s, p);
        #endif
    #else
        #ifdef use_pc
            Q_solver.solve(s,p_c);
            p += p_c;
        #else
            Q_solver.solve(s,p);
        #endif
    #endif
#else
    #ifdef use_pc
        p_c.v = Q_solver.solve(s.v);
        p += p_c;
    #else
        p.v = Q_solver.solve(s.v);
    #endif
#endif
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
#ifdef GPU_LIN_ALG
#ifdef testing_gpu
    solve(A, b.xv, v_n.xv);
    solve(A, b.yv, v_n.yv);
    solve(A, b.zv, v_n.zv);
#else
    A_solver.solve(b.xv, v_n.xv);
    A_solver.solve(b.yv, v_n.yv);
    A_solver.solve(b.zv, v_n.zv);
#endif
#else
    v_n.xv.v = A_solver.solve(b.xv.v);
    v_n.yv.v = A_solver.solve(b.yv.v);
    v_n.zv.v = A_solver.solve(b.zv.v);
#endif
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
        make_s(s, Re, dt, v_n, v_n1, p);
        //make_s_first(s, Re, dt, v_n);
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n" << std::flush;

        std::cout << "\tSolving for p" << std::flush;
        start = std::chrono::high_resolution_clock::now();
        //solve for p for the next timestep
#ifdef GPU_LIN_ALG
#ifdef testing_gpu
#ifdef use_pc
        solve(Q, s, p_c);
        p += p_c;
#else
        solve(Q, s, p);
#endif
#else
        #ifdef use_pc
            Q_solver.solve(s,p_c);
            p += p_c;
        #else
            Q_solver.solve(s,p);
        #endif
#endif
#else
        #ifdef use_pc
        p_c.v = Q_solver.solve(s.v);
        p += p_c;
    #else
        p.v = Q_solver.solve(s.v);
    #endif
#endif
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
#ifdef GPU_LIN_ALG
#ifdef testing_gpu
        solve(A, b.xv, v_n.xv);
#else
        A_solver.solve(b.xv, v_n.xv);
#endif
#else
        v_n.xv.v = A_solver.solve(b.xv.v);
#endif
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";


        std::cout << "\tSolving for vy" << std::flush;
        start = std::chrono::high_resolution_clock::now();
#ifdef GPU_LIN_ALG
#ifdef testing_gpu
        solve(A, b.yv, v_n.yv);
#else
        A_solver.solve(b.yv, v_n.yv);
#endif
#else
        v_n.yv.v = A_solver.solve(b.yv.v);
#endif
        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        std::cout << "\t" << dur.count() << "s\n";

        std::cout << "\tSolving for vz" << std::flush;
        start = std::chrono::high_resolution_clock::now();
#ifdef GPU_LIN_ALG
#ifdef testing_gpu
        solve(A, b.zv, v_n.zv);
#else
        A_solver.solve(b.zv, v_n.zv);
#endif
#else
        v_n.zv.v = A_solver.solve(b.zv.v);
#endif
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

    //write_vec(v_n, "../DEBUG/plotting_raw_data/end.txt");






}



#endif //CODE_CREATE_FLOW_HPP
