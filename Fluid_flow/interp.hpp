//
// Created by jacob on 22/10/21.
//

#ifndef CODE_INTERP_HPP
#define CODE_INTERP_HPP


#define STORE_SOLVERS
//#define STORE_MATS
#define LINEAR_INTERP

//#define LOG_INTERP_TIMES

namespace is {
#ifdef LINEAR_INTERP
    constexpr unsigned no_points = 8;
#else
    constexpr unsigned no_points = 24;//20;
#endif
}

typedef Eigen::FullPivLU<Eigen::Matrix<double, is::no_points, is::no_points> > eig_interp_solver;


template <typename T>
void get_interp_inds(T &interp_indices, const grid &g, unsigned index) {
    interp_indices[0] = index;

    for (const auto x : {1,-1}) {
        for (const auto y : {1,-1}) {
            for (const auto z : {1, -1}) {
                const bool is_good = g.can_move(index, x,0,0) &&g.can_move(index, 0,y,0) && g.can_move(index, 0,0,z) &&
                        g.can_move(index, x,y,0) && g.can_move(index, x,0,z) && g.can_move(index, 0,y,z)
                        && g.can_move(index, x,y,z);
                if (is_good) {
                    interp_indices[1] = g.get_move_ind(index, x, 0, 0);
                    interp_indices[2] = g.get_move_ind(index, 0, y, 0);
                    interp_indices[3] = g.get_move_ind(index, 0, 0, z);
                    interp_indices[4] = g.get_move_ind(index, x, y, 0);
                    interp_indices[5] = g.get_move_ind(index, x, 0, z);
                    interp_indices[6] = g.get_move_ind(index, 0, y, z);
                    interp_indices[7] = g.get_move_ind(index, x, y, z);
                    return;
                }

            }
        }
    }
}


#ifdef STORE_SOLVERS
struct interp_solvers {
                        template <class T>
                                class no_init_alloc
                                        : public std::allocator<T>
                                        {
                                        public:
                                            using std::allocator<T>::allocator;

                                            template <class U, class... Args> void construct(U*, Args&&...) {}
                                        };

                        std::vector<eig_interp_solver, no_init_alloc<eig_interp_solver>  >  solvers;
#ifdef STORE_MATS
                        std::vector<Eigen::Matrix<double, no_points, no_points>, no_init_alloc<Eigen::Matrix<double, no_points, no_points>> > mats;
#endif

                        explicit interp_solvers(unsigned no_grid_points) {
                            solvers.resize(no_grid_points);
#ifdef STORE_MATS
                            mats.resize(no_grid_points);
#endif
                        }


                        explicit interp_solvers(const grid &g) : interp_solvers(g.size()) {
                            //finding points to use for interpolation
                            // - assumes the offsets are smaller than the step size
        #pragma omp parallel for
        for (unsigned index = 0; index < g.size(); index++) {
            unsigned interp_indices[is::no_points];
            get_interp_inds(interp_indices, g, index);

            //finding the constants in y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
            //setting the matrix
#ifdef STORE_MATS

#else
Eigen::Matrix<double, is::no_points, is::no_points> mat;
#endif

for (unsigned j = 0; j < is::no_points; j++) {

    const auto x = g.x[interp_indices[j]];
    const auto y = g.y[interp_indices[j]];
    const auto z = g.z[interp_indices[j]];

#ifdef STORE_MATS
    solvers.mats[index](j, 0) = 1;
    solvers.mats[index](j, 1) = x;
    solvers.mats[index](j, 2) = y;
    solvers.mats[index](j, 3) = z;
    solvers.mats[index](j, 4) = x*y;
    solvers.mats[index](j, 5) = x*z;
    solvers.mats[index](j, 6) = y*z;
    solvers.mats[index](j, 7) = x*y*z;

    solvers.mats[index](j, 8) = x*x;
    solvers.mats[index](j, 9) = y*y;
    solvers.mats[index](j, 10) = z*z;
    solvers.mats[index](j, 11) = x*x*y;
    solvers.mats[index](j, 12) = x*x*z;
    solvers.mats[index](j, 13) = y*y*x;
    solvers.mats[index](j, 14) = y*y*z;
    solvers.mats[index](j, 15) = x*z*z;
    solvers.mats[index](j, 16) = y*z*z;
    solvers.mats[index](j, 17) = x*x*y*z;
    solvers.mats[index](j, 18) = x*y*y*z;
    solvers.mats[index](j, 19) = x*y*z*z;
    solvers.mats[index](j, 20) = x*x* y*y* z;
    solvers.mats[index](j, 21) = x*x *y *z*z;
    solvers.mats[index](j, 22) = x *y*y* z*z;
    solvers.mats[index](j,23) = x*x* y*y* z*z;
#else

#ifdef LINEAR_INTERP
    mat(j, 0) = 1;
    mat(j, 1) = x;
    mat(j, 2) = y;
    mat(j, 3) = z;
    mat(j, 4) = x*y;
    mat(j, 5) = x*z;
    mat(j, 6) = y*z;
    mat(j, 7) = x*y*z;
#else
    mat(j, 0) = 1;
    mat(j, 1) = x;
    mat(j, 2) = y;
    mat(j, 3) = z;
    mat(j, 4) = x*y;
    mat(j, 5) = x*z;
    mat(j, 6) = y*z;
    mat(j, 7) = x*y*z;

    mat(j, 8) = x*x;
    mat(j, 9) = y*y;
    mat(j, 10) = z*z;
    mat(j, 11) = x*x*y;
    mat(j, 12) = x*x*z;
    mat(j, 13) = y*y*x;
    mat(j, 14) = y*y*z;
    mat(j, 15) = x*z*z;
    mat(j, 16) = y*z*z;
    mat(j, 17) = x*x*y*z;
    mat(j, 18) = x*y*y*z;
    mat(j, 19) = x*y*z*z;
    mat(j, 20) = x*x* y*y* z;
    mat(j, 21) = x*x *y *z*z;
    mat(j, 22) = x *y*y* z*z;
    mat(j,23) = x*x* y*y* z*z;
#endif
#endif
}

#ifdef STORE_SOLVERS
#ifdef STORE_MATS
solvers.solvers[index] = Eigen::FullPivHouseholderQR<Eigen::Matrix<double, is::no_points, is::no_points> >(solvers.mats[index]);
#else
solvers[index] = eig_interp_solver(mat);
#endif
#else
const Eigen::FullPivHouseholderQR<Eigen::Matrix<double, is::no_points, is::no_points> > solver(mat);
#endif

        }


                        }


                    };
#endif




template <typename T>
double update_buffer(const T& a, const double x, const double y, const double z) {
#ifdef LINEAR_INTERP
    return a(0) + a(1)*x + a(2)*y + a(3)*z + a(4)*x*y + a(5)*x*z + a(6)*y*z + a(7)*x*y*z;
#else
    return a(0) + a(1)*x + a(2)*y + a(3)*z + a(4)*x*y + a(5)*x*z + a(6)*y*z + a(7)*x*y*z +
    a(8)*x*x + a(9)*y*y + a(10)*z*z +  a(11)*x*x*y + a(12)*x*x*z + a(13)*y*y*x +
    a(14)*y*y*z + a(15)*x*z*z + a(16)*y*z*z + a(17)*x*x*y*z + a(18)*x*y*y*z + a(19)*x*y*z*z +
    a(20)*x*x* y*y* z + a(21)*x*x *y *z*z + a(22)*x *y*y* z*z + a(23)*x*x* y*y* z*z;
#endif
}


#ifdef STORE_SOLVERS
void interpolate_vectors( big_vec_v &v_n, big_vec_v &v_n1, big_vec_d &p, const vec3& vel, const vec3& c_o_m, const vec3& omega, const double dt, interp_solvers & isolver, const grid &init_grid, const vec3& init_com) {
#else
        void interpolate_vectors( big_vec_v &v_n, big_vec_v &v_n1, big_vec_d &p, const vec3& vel, const vec3& c_o_m, const vec3& omega, const double dt) {
#endif

            const auto rot_angle_vec = omega*dt;
            const auto trans_vec = vel*dt;

            //variable to hold the new interpolated values
            Eigen::Matrix<double, Eigen::Dynamic, 1> vn_buff_x, vn_buff_y, vn_buff_z,
            vn1_buff_x, vn1_buff_y, vn1_buff_z, p_buff;
            vn_buff_x.resize(static_cast<long>(v_n.g->size()));
            vn_buff_y.resize(static_cast<long>(v_n.g->size()));
            vn_buff_z.resize(static_cast<long>(v_n.g->size()));

            vn1_buff_x.resize(static_cast<long>(v_n.g->size()));
            vn1_buff_y.resize(static_cast<long>(v_n.g->size()));
            vn1_buff_z.resize(static_cast<long>(v_n.g->size()));

            p_buff.resize(static_cast<long>(v_n.g->size()));





            //finding points to use for interpolation
            // - assumes the offsets are smaller than the step size

#ifdef LOG_INTERP_TIMES
            double time_finding_indices = 0;
            double time_filling_matrices = 0;
            double time_solving = 0;
            double time_moving = 0;
#ifdef STORE_SOLVERS

#else
            double time_solver = 0;
#endif
#endif
    #pragma omp parallel for
            for (unsigned index = 0; index < v_n.g->size(); index++) {
                unsigned interp_indices[is::no_points];

#ifdef STORE_SOLVERS
                const auto rot = rotate(init_com, rot_angle_vec, vec3(init_grid.x[index], init_grid.y[index], init_grid.z[index]), rot_angle_vec.length() );
#else
                const auto &g = *(v_n.g);
                const auto rot = rotate(c_o_m, rot_angle_vec, vec3(g.x[index], g.y[index], g.z[index]), rot_angle_vec.length() );
#endif

                const auto x = rot.x()+trans_vec.x();
                const auto y = rot.y()+trans_vec.y();
                const auto z = rot.z()+trans_vec.z();

                vec3 move_vec = {0,0,0};
#ifdef STORE_SOLVERS
                if (x < init_grid.x[index]) {
                    if (v_n.g->can_move(index, move_vec + vec3(-1,0,0)  )) {
                        move_vec += vec3(-1,0,0);
                    }
                }
                if (y < init_grid.y[index]) {
                    if (v_n.g->can_move(index, move_vec + vec3(0,-1,0)  )) {
                        move_vec +=vec3(0,-1,0);
                    }
                }
                if (z < init_grid.z[index]) {
                    if (v_n.g->can_move(index, move_vec + vec3(0,0,-1)  )) {
                        move_vec += vec3(0,0,-1);;
                    }
                }
#else
                if (x < g.x[index]) {
                    if (v_n.g->can_move(index, move_vec + vec3(-1,0,0)  )) {
                        move_vec += vec3(-1,0,0);
                    }
                }
                if (y < g.y[index]) {
                    if (v_n.g->can_move(index, move_vec + vec3(0,-1,0)  )) {
                        move_vec +=vec3(0,-1,0);
                    }
                }
                if (z < g.z[index]) {
                    if (v_n.g->can_move(index, move_vec + vec3(0,0,-1)  )) {
                        move_vec += vec3(0,0,-1);;
                    }
                }
#endif
                const auto i_index = v_n.g->get_move_ind(index, move_vec);

#ifdef LOG_INTERP_TIMES
                const auto start_indices = std::chrono::high_resolution_clock::now();
#endif

                get_interp_inds(interp_indices, *v_n.g, i_index);

#ifdef LOG_INTERP_TIMES
                const auto end_indices = std::chrono::high_resolution_clock::now();
                time_finding_indices += static_cast<std::chrono::duration<double>>(end_indices - start_indices).count();
#endif





                //finding the constants in y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
                //setting the matrix
#ifdef LOG_INTERP_TIMES
                const auto start_filling = std::chrono::high_resolution_clock::now();
#endif

#ifdef STORE_MATS

#else
                Eigen::Matrix<double, is::no_points, is::no_points> mat;
#endif

                Eigen::Matrix<double, is::no_points, 1> vn_vec_x, vn_vec_y, vn_vec_z,
                vn1_vec_x, vn1_vec_y, vn1_vec_z, p_vec;




                for (unsigned j = 0; j < is::no_points; j++) {
                    vn_vec_x[j] = v_n.xv.v[interp_indices[j]];
                    vn_vec_y[j] = v_n.yv.v[interp_indices[j]];
                    vn_vec_z[j] = v_n.zv.v[interp_indices[j]];

                    vn1_vec_x[j] = v_n1.xv.v[interp_indices[j]];
                    vn1_vec_y[j] = v_n1.yv.v[interp_indices[j]];
                    vn1_vec_z[j] = v_n1.zv.v[interp_indices[j]];

                    p_vec[j] = p.v[interp_indices[j]];



#ifdef STORE_SOLVERS

#else
                    const auto x = g.x[interp_indices[j]];
                    const auto y = g.y[interp_indices[j]];
                    const auto z = g.z[interp_indices[j]];

                    mat(j, 0) = 1;
                    mat(j, 1) = x;
                    mat(j, 2) = y;
                    mat(j, 3) = z;
                    mat(j, 4) = x*y;
                    mat(j, 5) = x*z;
                    mat(j, 6) = y*z;
                    mat(j, 7) = x*y*z;

                    /*mat(j, 8) = x*x;
                    mat(j, 9) = y*y;
                    mat(j, 10) = z*z;
                    mat(j, 11) = x*x*y;
                    mat(j, 12) = x*x*z;
                    mat(j, 13) = y*y*x;
                    mat(j, 14) = y*y*z;
                    mat(j, 15) = x*z*z;
                    mat(j, 16) = y*z*z;
                    mat(j, 17) = x*x*y*z;
                    mat(j, 18) = x*y*y*z;
                    mat(j, 19) = x*y*z*z;
                    mat(j, 20) = x*x* y*y* z;
                    mat(j, 21) = x*x *y *z*z;
                    mat(j, 22) = x *y*y* z*z;
                    mat(j,23) = x*x* y*y* z*z;*/
#endif
                }
#ifdef LOG_INTERP_TIMES
                const auto end_filling = std::chrono::high_resolution_clock::now();
                time_filling_matrices += static_cast<std::chrono::duration<double>>(end_filling - start_filling).count();
#endif


#ifdef STORE_SOLVERS
#else
                const auto start_solver = std::chrono::high_resolution_clock::now();
#endif

#ifdef STORE_SOLVERS

#else
                const Eigen::FullPivHouseholderQR<Eigen::Matrix<double, is::no_points, is::no_points> > solver(mat);
#endif

#ifdef STORE_SOLVERS
#else
        #ifdef LOG_INTERP_TIMES
                const auto end_solver = std::chrono::high_resolution_clock::now();
                time_solver += static_cast<std::chrono::duration<double>>(end_solver - start_solver).count();
        #endif
#endif

                const auto start_solving = std::chrono::high_resolution_clock::now();
                //mat and vec now set, just need to solve for the coefficients
                //solver.compute(mat);
#ifdef STORE_SOLVERS
const decltype(vn_vec_x) a_vn_x = isolver.solvers[i_index].solve(vn_vec_x);
const decltype(vn_vec_y) a_vn_y = isolver.solvers[i_index].solve(vn_vec_y);
const decltype(vn_vec_z) a_vn_z = isolver.solvers[i_index].solve(vn_vec_z);
const decltype(vn1_vec_x) a_vn1_x = isolver.solvers[i_index].solve(vn1_vec_x);
const decltype(vn1_vec_y) a_vn1_y = isolver.solvers[i_index].solve(vn1_vec_y);
const decltype(vn1_vec_z) a_vn1_z = isolver.solvers[i_index].solve(vn1_vec_z);
const decltype(p_vec) a_p = isolver.solvers[i_index].solve(p_vec);
#else
const decltype(vn_vec_x) a_vn_x = solver.solve(vn_vec_x);
const decltype(vn_vec_y) a_vn_y = solver.solve(vn_vec_y);
const decltype(vn_vec_z) a_vn_z = solver.solve(vn_vec_z);
const decltype(vn1_vec_x) a_vn1_x = solver.solve(vn1_vec_x);
const decltype(vn1_vec_y) a_vn1_y = solver.solve(vn1_vec_y);
const decltype(vn1_vec_z) a_vn1_z = solver.solve(vn1_vec_z);
const decltype(p_vec) a_p = solver.solve(p_vec);
#endif


#ifndef NDEBUG
if ( std::abs(x-v_n.g->x[index]) > v_n.g->dx) {
    std::cerr << "movement in the x-direction is larger than step size. Interpolation of values might be inaccurate\n";
}
if (std::abs(y-v_n.g->y[index]) > v_n.g->dy) {
    std::cerr << "movement in the y-direction is larger than step size. Interpolation of values might be inaccurate\n";
}
if (std::abs(z-v_n.g->z[index]) > v_n.g->dz) {
    std::cerr << "movement in the z-direction is larger than step size. Interpolation of values might be inaccurate\n";
}
#endif

vn_buff_x[index] = update_buffer(a_vn_x, x,y,z);
vn_buff_y[index] = update_buffer(a_vn_y, x,y,z);
vn_buff_z[index] = update_buffer(a_vn_z, x,y,z);
vn1_buff_x[index] = update_buffer(a_vn1_x, x,y,z);
vn1_buff_y[index] = update_buffer(a_vn1_y, x,y,z);
vn1_buff_z[index] = update_buffer(a_vn1_z, x,y,z);
p_buff[index] = update_buffer(a_p, x,y,z);

#ifdef LOG_INTERP_TIMES
const auto end_solving = std::chrono::high_resolution_clock::now();
time_solving += static_cast<std::chrono::duration<double>>(end_solving - start_solving).count();
#endif

            }
#ifdef LOG_INTERP_TIMES
            const auto start_moving = std::chrono::high_resolution_clock::now();
#endif
            //v = std::move(buffer);
            v_n.xv.v = std::move(vn_buff_x);
            v_n.yv.v = std::move(vn_buff_y);
            v_n.zv.v = std::move(vn_buff_z);

            v_n1.xv.v = std::move(vn1_buff_x);
            v_n1.yv.v = std::move(vn1_buff_y);
            v_n1.zv.v = std::move(vn1_buff_z);

            p.v = std::move(p_buff);
#ifdef LOG_INTERP_TIMES
            const auto end_moving = std::chrono::high_resolution_clock::now();
            time_moving += static_cast<std::chrono::duration<double>>(end_moving - start_moving).count();

            std::cerr << "\n\tTime spent finding indices : " << time_finding_indices << "\n";
            std::cerr << "\tTime spent filling matrices : " << time_filling_matrices << "\n";
#ifdef STORE_SOLVERS
#else
            std::cerr << "\tTime spent initialising solver : " << time_solver << "\n";
#endif
            std::cerr << "\tTime spent solving system : " << time_solving << "\n";
            std::cerr << "\tTime spent moving : " << time_moving << "\n";
#endif
        }








#endif //CODE_INTERP_HPP
