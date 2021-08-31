//
// Created by jacob on 31/8/21.
//

#ifndef CODE_UPDATE_MESH_HPP
#define CODE_UPDATE_MESH_HPP

#include "boundary_conditions.hpp"
#include "../Rigid_body/body.hpp"

template <unsigned N, unsigned M, unsigned P>
void update_mesh(const boundary_conditions<N,M,P> &bc, body *b, double dt) {
    const auto dx = bc.p_bc.dx;
    const auto dy = bc.p_bc.dy;
    const auto dz = bc.p_bc.dz;


    std::vector<vec3> forces, points;
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
                    if (bc.norms.contains(i,j,k) && bc.norms.normal(i,j,k) != vec3(0)) {
                        points[forces_counter] = vec3(i/dx+dx/2, j/dy+dy/2, k/dz+dz/2); //forces applied at the center of grid points
                        //taking the force as pointing against the normal, not sure if this is right
                        forces[forces_counter++] = bc.p_bc(i,j,k)*dx*dy*dz * - bc.norms.normal(i,j,k); //P=F/A  =>  F=PA
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

#endif //CODE_UPDATE_MESH_HPP
