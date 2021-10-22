//
// Created by jacob on 31/8/21.
//

#ifndef CODE_UPDATE_MESH_HPP
#define CODE_UPDATE_MESH_HPP

#include "boundary_conditions.hpp"
#include "../Rigid_body/body.hpp"
#include "update_vecs.hpp"
#include "interp.hpp"

#ifdef FLUID_MOVES_MESH
bool fluid_moves(const double t) {
    return t>0.01;
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

    forces.resize( bc.size() );
    points.resize( forces.size() );

#ifndef NDEBUG
        const auto forces_is = forces.size();
#endif
        unsigned forces_counter = 0;

        //could be done more efficiently using a range base for loop through bc.m_points
        for (unsigned i = 0; i < p.size(); i++) {
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

#ifdef STORE_SOLVERS
    interpolate_vectors(v_n, v_n1, p, b->model.v, old_c_o_m, b->model.w, dt, isolver, init_grid, init_com);
#else
interpolate_vectors(v_n, v_n1, p, b->model.v, old_c_o_m, b->model.w, dt);
#endif

    v_n.g->update_pos(b->model.v, old_c_o_m, b->model.w, dt);

    bc.update(b->model.w, b->model.v, old_c_o_m, dt);


    if (!update_pressure_BC<false>(bc, p)) {
        throw std::runtime_error("enforcing pressure boundary condition failed");
    }


    //can't think of a better way to make sure that the extrapolation does not affect points that need to have BC enforced
    if (!enforce_velocity_BC<false>(bc, v_n)) {
        throw std::runtime_error("enforcing velocity boundary condition failed");
    }


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

#endif //CODE_UPDATE_MESH_HPP
