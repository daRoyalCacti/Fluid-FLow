//
// Created by jacob on 31/8/21.
//

#ifndef CODE_UPDATE_MESH_HPP
#define CODE_UPDATE_MESH_HPP

#include "boundary_conditions.hpp"
#include "../Rigid_body/body.hpp"

bool fluid_moves(const double t) {
    return false;
}

vec3 global_forces(const double t) {
    return vec3(0);
}

//updating mesh also requires updating the vectors
// - have to extrapolate p, v_n but also v_n1 because some equations require it
template <unsigned N, unsigned M, unsigned P>
void update_mesh(boundary_conditions<N,M,P> &bc, body *b, const big_vec<N,M,P, vec3> &v_n, const big_vec<N,M,P, vec3> &v_n1, const big_vec<N,M,P, double> &p, const double dt, const double t) {
    const auto dx = bc.p_bc.dx;
    const auto dy = bc.p_bc.dy;
    const auto dz = bc.p_bc.dz;


    std::vector<vec3> forces, points;
    if (fluid_moves(t)) {
        forces.resize( bc.norms.size() - bc.no_wall_points() - bc.no_inside_mesh );
        points.resize( forces.size() );

#ifndef NDEBUG
        const auto forces_is = forces.size();
#endif
        unsigned forces_counter = 0;

        for (unsigned i = 0; i <= N; i++) {
            for (unsigned j = 0; j <= M; j++) {
                for (unsigned k = 0; k <= P; k++) {
                    if (i > 0 && i <N && j > 0 && j< M && k>0 && k<P) {  //if off the boundary
                        if (bc.norms.contains(i,j,k) && bc.norms.normal(i,j,k) != vec3(0)) {    //if boundary point outside of mesh
                            points[forces_counter] = vec3(i/dx+dx/2, j/dy+dy/2, k/dz+dz/2); //forces applied at the center of grid points
                            //taking the force as pointing against the normal, not sure if this is right
                            forces[forces_counter++] = bc.p_bc(i,j,k)*dx*dy*dz * - bc.norms.normal(i,j,k)  +
                                    global_forces(t); //P=F/A  =>  F=PA
                        }
                    }
                }
            }
        }

#ifndef NDEBUG
        if (forces_counter != forces_is) {
            std::cerr << "forces is the wrong size for the number of boundary points\n";
        }
#endif


        b->update_pos(forces, points, dt);
    }
    bc.extrapolate(v_n);
    bc.extrapolate(v_n1);
    bc.extrapolate(p);
    bc.update_mesh_boundary();
}


template <unsigned N, unsigned M, unsigned P>
void update(boundary_conditions<N,M,P> &bc, body *b, const double dt, const double t) {
    update_mesh(bc, b, dt, t);
}

#endif //CODE_UPDATE_MESH_HPP
