//
// Created by jacob on 31/8/21.
//

#ifndef CODE_UPDATE_MESH_HPP
#define CODE_UPDATE_MESH_HPP

#include "boundary_conditions.hpp"
#include "../Rigid_body/body.hpp"
#include "update_vecs.hpp"

#ifdef FLUID_MOVES_MESH
bool fluid_moves(const double t) {
    return t>0.001; //0.01
}
#else
bool fluid_moves(const double t) {
    return false;
}
#endif

vec3 global_forces(const double t) {
    //return {0.00001*sin(t), 0, 0};
    return vec3(0);
}

#define STORE_SOLVERS
//#define STORE_MATS
#define LINEAR_INTERP

namespace is {
#ifdef LINEAR_INTERP
    constexpr unsigned no_points = 8;
#else
    constexpr unsigned no_points = 24;//20;
#endif
}

//https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
// - checked in mathematica, matrix is not definite, negative semi-definite or positive semi-definite
//   * means LLT and LDLT are out
// - matrix is invertible (checked in mathematica)
// - really want +++ interms of accuracy. Means its a choice between
//   *(+)CompleteOrthogonalDecomposition (took 9.08736 - although first real timestep too only 5.2787)
//   *(-)FullPivHouseholderQR (took 9.98239)
//   *(+)ColPivHouseholderQR (took 4.78679)
//   *(-)FullPivLU (took 2.73205)
//   *(-)BDCSVD (gave infs immediately)
//   *(-)JacobiSVD (gave infs immediately)
//          /\ timings for second real timestep with linear interpolation
typedef Eigen::FullPivLU<Eigen::Matrix<double, is::no_points, is::no_points> > eig_interp_solver;


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


        double time_finding_indices = 0;
        double time_filling_matrices = 0;
        double time_solver = 0;

    #pragma omp parallel for
        for (unsigned index = 0; index < g.size(); index++) {
            unsigned interp_indices[is::no_points];
            size_t counter = 0;

            const auto start_indices = std::chrono::high_resolution_clock::now();
            //include the point itself
            interp_indices[counter++] = index;
            int i = 1;
            while (true) {

                //getting corner values first
                for (const auto& horiz : {-i,i}) {
                    for (const auto& vert : {-i,i}) {
                        for (const auto& in : {-i,i}) {
                            if (g.can_move(index, horiz, vert, in)) {
                                interp_indices[counter++] = g.get_move_ind(index, horiz, vert, in);
                                if (counter == is::no_points) {
                                    goto got_indices;
                                }
                            }
                        }
                    }
                }

                //then getting edge values
                for (const auto& horiz : {-i,i}) {
                    for (const auto& vert : {-i,i}) {
                        if (g.can_move(index, horiz, vert, 0)) {
                            interp_indices[counter++] = g.get_move_ind(index, horiz, vert, 0);
                            if (counter == is::no_points) {
                                goto got_indices;
                            }
                        }
                    }

                    for (const auto& in : {-i,i}) {
                        if (g.can_move(index, horiz, 0, in)) {
                            interp_indices[counter++] = g.get_move_ind(index, horiz, 0, in);
                            if (counter == is::no_points) {
                                goto got_indices;
                            }
                        }
                    }
                }


                //then points along major axes
                for (const auto& horiz : {-i,i}) {
                    if (g.can_move(index, horiz, 0, 0)) {
                        interp_indices[counter++] = g.get_move_ind(index, horiz, 0, 0);
                        if (counter == is::no_points) {
                            goto got_indices;
                        }
                    }
                }
                for (const auto& vert : {-i,i}) {
                    if (g.can_move(index, 0, vert, 0)) {
                        interp_indices[counter++] = g.get_move_ind(index, 0, vert, 0);
                        if (counter == is::no_points) {
                            goto got_indices;
                        }
                    }
                }
                for (const auto& in : {-i,i}) {
                    if (g.can_move(index, 0, 0, in)) {
                        interp_indices[counter++] = g.get_move_ind(index, 0, 0, in);
                        if (counter == is::no_points) {
                            goto got_indices;
                        }
                    }
                }

                //getting the rest of the edge values
                for (const auto& in : {-i,i}) {
                    for (const auto& vert : {-i,i}) {
                        if (g.can_move(index, 0, vert, in)) {
                            interp_indices[counter++] = g.get_move_ind(index, 0, vert, in);
                            if (counter == is::no_points) {
                                goto got_indices;
                            }
                        }
                    }
                }




                i++;
            }
            got_indices:
            const auto end_indices = std::chrono::high_resolution_clock::now();
            time_finding_indices += static_cast<std::chrono::duration<double>>(end_indices - start_indices).count();


            //finding the constants in y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
            //setting the matrix
            const auto start_filling = std::chrono::high_resolution_clock::now();
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

            const auto end_filling = std::chrono::high_resolution_clock::now();
            time_filling_matrices += static_cast<std::chrono::duration<double>>(end_filling - start_filling).count();



            const auto start_solver = std::chrono::high_resolution_clock::now();
            //const Eigen::LDLT<Eigen::Matrix<double, no_points, no_points> > solver(mat);  //bad
            //const Eigen::FullPivLU<Eigen::Matrix<double, no_points, no_points> > solver(mat); //maybe bad
            //const Eigen::FullPivHouseholderQR<Eigen::Matrix<double, no_points, no_points> > solver(mat);
#ifdef STORE_SOLVERS
#ifdef STORE_MATS
solvers.solvers[index] = Eigen::FullPivHouseholderQR<Eigen::Matrix<double, is::no_points, is::no_points> >(solvers.mats[index]);
#else
solvers[index] = eig_interp_solver(mat);
#endif
#else
const Eigen::FullPivHouseholderQR<Eigen::Matrix<double, is::no_points, is::no_points> > solver(mat);
#endif

const auto end_solver = std::chrono::high_resolution_clock::now();
time_solver += static_cast<std::chrono::duration<double>>(end_solver - start_solver).count();
        }
        std::cerr << "\n\tTime spent finding indices : " << time_finding_indices << "\n";
        std::cerr << "\tTime spent filling matrices : " << time_filling_matrices << "\n";
        std::cerr << "\tTime spent initialising solver : " << time_solver << "\n";

    }


};
#endif


//defined below update_mesh
#ifdef STORE_SOLVERS
void interpolate_vectors( big_vec_v &v_n, big_vec_v &v_n1, big_vec_d &p, const vec3& vel, const vec3& c_o_m, const vec3& omega, double dt, interp_solvers & isolver, const grid &init_grid, const vec3& init_com);
#else
void interpolate_vectors( big_vec_v &v_n, big_vec_v &v_n1, big_vec_d &p, const vec3& vel, const vec3& c_o_m, const vec3& omega, double dt);
#endif

//updating mesh also requires updating the vectors
// - have to extrapolate p, v_n but also v_n1 because some equations require it
//counter just for debug
#ifdef STORE_SOLVERS
void update_mesh(boundary_conditions &bc, body *b, big_vec_v &v_n, big_vec_v &v_n1, big_vec_d &p, const double dt, const double t, interp_solvers & isolver, const grid &init_grid, const vec3& init_com, const unsigned counter = 0) {
#else
void update_mesh(boundary_conditions &bc, body *b, big_vec_v &v_n, big_vec_v &v_n1, big_vec_d &p, const double dt, const double t, const unsigned counter = 0) {
#endif
    if (!enforce_velocity_BC<false>(bc, v_n)) {
        throw std::runtime_error("enforcing velocity boundary condition failed");
    }
    if (!update_pressure_BC<false>(bc, p) ) {
        throw std::runtime_error("enforcing pressure boundary condition failed");
    }

    std::vector<vec3> forces, points;
    const auto& g = *v_n.g;

    //forces.resize( bc.norms.size() - bc.no_wall_points() - bc.no_inside_mesh );
    forces.resize( bc.size() );
    points.resize( forces.size() );

#ifndef NDEBUG
        const auto forces_is = forces.size();
#endif
        unsigned forces_counter = 0;

        //could be done more efficiently using a range base for loop through bc.m_points
        for (unsigned i = 0; i < p.size(); i++) {
            //std::cerr << i << "\n";
            if (g.off_walls(i)) {  //if off the boundary
                if (g.is_boundary(i)) {    //if boundary point
                    points[forces_counter] = bc.m_points.get_point(i);    //forces applied at the mesh boundary
                    //taking the force as pointing against the normal, not sure if this is right
                    if (fluid_moves(t)) {
                        const auto dx = g.dx;
                        const auto dy = g.dy;
                        const auto dz = g.dz;
                        forces[forces_counter++] = p(i)*dx*dy*dz * - bc.norms.normal(i)  +
                                global_forces(t); //P=F/A  =>  F=PA
                    } else {
                        forces[forces_counter++] = global_forces(t);
                    }
                }

            }
        }

#ifndef NDEBUG
        if (forces_counter != forces_is) {
            std::cerr << "forces is the wrong size for the number of boundary points\n";
            std::cerr << "size of forces : " << forces_is << ". Number of boundary points : " << forces_counter << "\n";
        }
#endif

    const auto old_c_o_m = b->model.pos_cm;
    b->update_pos(forces, points, dt);

    //extrapolate must be called before update_mesh_boundary because this requires to old values for the normals
    /*bc.extrapolate(v_n);
    bc.extrapolate(v_n1);
    bc.extrapolate(p);*/

    //bc.enforce_velocity_BC(v_n1);
    /*bc.update_mesh_boundary();*/

    /*std::cout << "\tUpdating v_n\n";
    v_n.move(b->model.v.x()*dt, b->model.v.y()*dt, b->model.v.z()*dt);
    std::cout << "\tUpdating v_n1\n";
    v_n1.move(b->model.v.x()*dt, b->model.v.y()*dt, b->model.v.z()*dt);
    std::cout << "\tUpdating p\n";
    p.move(b->model.v.x()*dt, b->model.v.y()*dt, b->model.v.z()*dt);
    std::cout << "\tUpdating g\n";
    v_n.g->move(b->model.v.x()*dt, b->model.v.y()*dt, b->model.v.z()*dt);   //only need to update g for one of the vectors since they all point to the same place*/
    //std::cerr << vec3(b->model.v.x()*dt, b->model.v.y()*dt, b->model.v.z()*dt) << "\n";

    /*const auto x_off = b->model.v.x()*dt;
    const auto y_off = b->model.v.y()*dt;
    const auto z_off = b->model.v.z()*dt;*/
    //const vec3& vel, const vec3& c_o_m, const vec3& omega, const double dt
#ifdef STORE_SOLVERS
interpolate_vectors(v_n, v_n1, p, b->model.v, old_c_o_m, b->model.w, dt, isolver, init_grid, init_com);
#else
interpolate_vectors(v_n, v_n1, p, b->model.v, old_c_o_m, b->model.w, dt);
#endif

    v_n.g->update_pos(b->model.v, old_c_o_m, b->model.w, dt);

    bc.update(b->model.w, b->model.v, old_c_o_m, dt);


    //std::cout << "\tenforcing boundary conditions\n";
    if (!update_pressure_BC<false>(bc, p)) {
        throw std::runtime_error("enforcing pressure boundary condition failed");
    }

    //updating pressure and wall velocity points
    // - non-wall velocity points and normals have already been updated
    //bc.update_velocity_wall_BC();
    //bc.update_pressure_BC(p);

    //update_wall_velocity(bc, v_n);


    //can't think of a better way to make sure that the extrapolation does not affect points that need to have BC enforced
    if (!enforce_velocity_BC<false>(bc, v_n)) {
        throw std::runtime_error("enforcing velocity boundary condition failed");
    }
    //bc.enforce_pressure_BC(p);

    std::cerr << "writing i files\n";
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

    const auto inds = v_n.g->get_middle_inds();
    write_vec(v_n,inds, ("../DEBUG/velocity_interpolated/" + file_name + ".txt").data());
    write_vec(p,inds, ("../DEBUG/pressure_interpolated/" + file_name + ".txt").data());
}

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



//const vec3& vel, const vec3& c_o_m, const vec3& omega, const double dt
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


    double time_finding_indices = 0;
    double time_filling_matrices = 0;
    double time_solving = 0;
    double time_moving = 0;
#ifdef STORE_SOLVERS

#else
    double time_solver = 0;
#endif

    #pragma omp parallel for
    for (unsigned index = 0; index < v_n.g->size(); index++) {
        unsigned interp_indices[is::no_points];
        size_t counter = 0;

        const auto start_indices = std::chrono::high_resolution_clock::now();
        //include the point itself
        interp_indices[counter++] = index;
        int i = 1;
        while (true) {

            //std::cerr << "\t" << counter << "\n";
            //getting corner values first
            for (const auto& horiz : {-i,i}) {
                for (const auto& vert : {-i,i}) {
                    for (const auto& in : {-i,i}) {
                        if (v_n.g->can_move(index, horiz, vert, in)) {
                            interp_indices[counter++] = v_n.g->get_move_ind(index, horiz, vert, in);
                            if (counter == is::no_points) {
                                goto got_indices;
                            }
                        }
                    }
                }
            }

            //then getting edge values
            for (const auto& horiz : {-i,i}) {
                for (const auto& vert : {-i,i}) {
                    if (v_n.g->can_move(index, horiz, vert, 0)) {
                        interp_indices[counter++] = v_n.g->get_move_ind(index, horiz, vert, 0);
                        if (counter == is::no_points) {
                            goto got_indices;
                        }
                    }
                }

                for (const auto& in : {-i,i}) {
                    if (v_n.g->can_move(index, horiz, 0, in)) {
                        interp_indices[counter++] = v_n.g->get_move_ind(index, horiz, 0, in);
                        if (counter == is::no_points) {
                            goto got_indices;
                        }
                    }
                }
            }


            //then points along major axes
            for (const auto& horiz : {-i,i}) {
                if (v_n.g->can_move(index, horiz, 0, 0)) {
                    interp_indices[counter++] = v_n.g->get_move_ind(index, horiz, 0, 0);
                    if (counter == is::no_points) {
                        goto got_indices;
                    }
                }
            }
            for (const auto& vert : {-i,i}) {
                if (v_n.g->can_move(index, 0, vert, 0)) {
                    interp_indices[counter++] = v_n.g->get_move_ind(index, 0, vert, 0);
                    if (counter == is::no_points) {
                        goto got_indices;
                    }
                }
            }
            for (const auto& in : {-i,i}) {
                if (v_n.g->can_move(index, 0, 0, in)) {
                    interp_indices[counter++] = v_n.g->get_move_ind(index, 0, 0, in);
                    if (counter == is::no_points) {
                        goto got_indices;
                    }
                }
            }

            //getting the rest of the edge values
            for (const auto& in : {-i,i}) {
                for (const auto& vert : {-i,i}) {
                    if (v_n.g->can_move(index, 0, vert, in)) {
                        interp_indices[counter++] = v_n.g->get_move_ind(index, 0, vert, in);
                        if (counter == is::no_points) {
                            goto got_indices;
                        }
                    }
                }
            }




            i++;
        }
        got_indices:
        const auto end_indices = std::chrono::high_resolution_clock::now();
        time_finding_indices += static_cast<std::chrono::duration<double>>(end_indices - start_indices).count();





        //finding the constants in y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
        //setting the matrix
        const auto start_filling = std::chrono::high_resolution_clock::now();
#ifdef STORE_MATS

#else
        Eigen::Matrix<double, is::no_points, is::no_points> mat;
#endif
        /*Eigen::BiCGSTAB<Eigen::Matrix<double, no_points, no_points> > solver;
        solver.setTolerance(1e-4);*/



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
        }

        const auto end_filling = std::chrono::high_resolution_clock::now();
        time_filling_matrices += static_cast<std::chrono::duration<double>>(end_filling - start_filling).count();


#ifdef STORE_SOLVERS
#else
        const auto start_solver = std::chrono::high_resolution_clock::now();
#endif
        //const Eigen::LDLT<Eigen::Matrix<double, no_points, no_points> > solver(mat);  //bad
        //const Eigen::FullPivLU<Eigen::Matrix<double, no_points, no_points> > solver(mat); //maybe bad
        //const Eigen::FullPivHouseholderQR<Eigen::Matrix<double, no_points, no_points> > solver(mat);
#ifdef STORE_SOLVERS

#else
        const Eigen::FullPivHouseholderQR<Eigen::Matrix<double, is::no_points, is::no_points> > solver(mat);
#endif

#ifdef STORE_SOLVERS
#else
        const auto end_solver = std::chrono::high_resolution_clock::now();
        time_solver += static_cast<std::chrono::duration<double>>(end_solver - start_solver).count();
#endif

        const auto start_solving = std::chrono::high_resolution_clock::now();
        //mat and vec now set, just need to solve for the coefficients
        //solver.compute(mat);
#ifdef STORE_SOLVERS
        const decltype(vn_vec_x) a_vn_x = isolver.solvers[index].solve(vn_vec_x);
        const decltype(vn_vec_y) a_vn_y = isolver.solvers[index].solve(vn_vec_y);
        const decltype(vn_vec_z) a_vn_z = isolver.solvers[index].solve(vn_vec_z);
        const decltype(vn1_vec_x) a_vn1_x = isolver.solvers[index].solve(vn1_vec_x);
        const decltype(vn1_vec_y) a_vn1_y = isolver.solvers[index].solve(vn1_vec_y);
        const decltype(vn1_vec_z) a_vn1_z = isolver.solvers[index].solve(vn1_vec_z);
        const decltype(p_vec) a_p = isolver.solvers[index].solve(p_vec);
#else
        const decltype(vn_vec_x) a_vn_x = solver.solve(vn_vec_x);
        const decltype(vn_vec_y) a_vn_y = solver.solve(vn_vec_y);
        const decltype(vn_vec_z) a_vn_z = solver.solve(vn_vec_z);
        const decltype(vn1_vec_x) a_vn1_x = solver.solve(vn1_vec_x);
        const decltype(vn1_vec_y) a_vn1_y = solver.solve(vn1_vec_y);
        const decltype(vn1_vec_z) a_vn1_z = solver.solve(vn1_vec_z);
        const decltype(p_vec) a_p = solver.solve(p_vec);
#endif

        /*const decltype(vn_vec_x) a_vn_x = mat.lu().solve(vn_vec_x);
        const decltype(vn_vec_y) a_vn_y = mat.lu().solve(vn_vec_y);
        const decltype(vn_vec_z) a_vn_z = mat.lu().solve(vn_vec_z);
        const decltype(vn1_vec_x) a_vn1_x = mat.lu().solve(vn_vec_x);
        const decltype(vn1_vec_y) a_vn1_y = mat.lu().solve(vn1_vec_y);
        const decltype(vn1_vec_z) a_vn1_z = mat.lu().solve(vn1_vec_z);
        const decltype(p_vec) a_p = mat.lu().solve(p_vec);*/


        /*const decltype(mat) mat_inv = mat.inverse();
        const decltype(vn_vec_x) a_vn_x = mat_inv*(vn_vec_x);
        const decltype(vn_vec_y) a_vn_y = mat_inv*(vn_vec_y);
        const decltype(vn_vec_z) a_vn_z = mat_inv*(vn_vec_z);
        const decltype(vn1_vec_x) a_vn1_x = mat_inv*(vn_vec_x);
        const decltype(vn1_vec_y) a_vn1_y = mat_inv*(vn1_vec_y);
        const decltype(vn1_vec_z) a_vn1_z = mat_inv*(vn1_vec_z);
        const decltype(p_vec) a_p = mat_inv*(p_vec);*/

        /*const auto rot_angle_vec = omega*dt;
        const auto trans_vec = vel*dt;

        for (unsigned i = 0; i < x.size(); i++) {
            const auto rot = rotate(c_o_m, rot_angle_vec, vec3(x[i], y[i], z[i]), rot_angle_vec.length() );
            x[i] = rot.x()+trans_vec.x();
            y[i] = rot.y()+trans_vec.y();
            z[i] = rot.z()+trans_vec.z();
        }*/

        /*const auto x = g.x[index] + x_off;
        const auto y = g.y[index] + y_off;
        const auto z = g.z[index] + z_off;*/
#ifdef STORE_SOLVERS
        const auto rot = rotate(init_com, rot_angle_vec, vec3(init_grid.x[index], init_grid.y[index], init_grid.z[index]), rot_angle_vec.length() );
#else
        const auto rot = rotate(c_o_m, rot_angle_vec, vec3(g.x[index], g.y[index], g.z[index]), rot_angle_vec.length() );
#endif
        const auto x = rot.x()+trans_vec.x();
        const auto y = rot.y()+trans_vec.y();
        const auto z = rot.z()+trans_vec.z();

#ifndef NDEBUG
        if ( std::abs(x-g.x[index]) > g.dx) {
            std::cerr << "movement in the x-direction is larger than step size. Interpolation of values might be inaccurate\n";
        }
        if (std::abs(y-g.y[index]) > g.dy) {
            std::cerr << "movement in the y-direction is larger than step size. Interpolation of values might be inaccurate\n";
        }
        if (std::abs(z-g.z[index]) > g.dz) {
            std::cerr << "movement in the z-direction is larger than step size. Interpolation of values might be inaccurate\n";
        }
#endif


        /*buffer[index] = a(0) + a(1)*x + a(2)*y + a(3)*z + a(4)*x*y + a(5)*x*z + a(6)*y*z + a(7)*x*y*z +
                a(8)*x*x + a(9)*y*y + a(10)*z*z +  a(11)*x*x*y + a(12)*x*x*z + a(13)*y*y*x +
                a(14)*y*y*z + a(15)*x*z*z + a(16)*y*z*z + a(17)*x*x*y*z + a(18)*x*y*y*z + a(19)*x*y*z*x ;*/
        vn_buff_x[index] = update_buffer(a_vn_x, x,y,z);
        vn_buff_y[index] = update_buffer(a_vn_y, x,y,z);
        vn_buff_z[index] = update_buffer(a_vn_z, x,y,z);
        vn1_buff_x[index] = update_buffer(a_vn1_x, x,y,z);
        vn1_buff_y[index] = update_buffer(a_vn1_y, x,y,z);
        vn1_buff_z[index] = update_buffer(a_vn1_z, x,y,z);
        p_buff[index] = update_buffer(a_p, x,y,z);

        const auto end_solving = std::chrono::high_resolution_clock::now();
        time_solving += static_cast<std::chrono::duration<double>>(end_solving - start_solving).count();

    }

    const auto start_moving = std::chrono::high_resolution_clock::now();
    //v = std::move(buffer);
    v_n.xv.v = std::move(vn_buff_x);
    v_n.yv.v = std::move(vn_buff_y);
    v_n.zv.v = std::move(vn_buff_z);

    v_n1.xv.v = std::move(vn1_buff_x);
    v_n1.yv.v = std::move(vn1_buff_y);
    v_n1.zv.v = std::move(vn1_buff_z);

    p.v = std::move(p_buff);

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
}



#endif //CODE_UPDATE_MESH_HPP
