//
// Created by jacob on 31/8/21.
//

#ifndef CODE_UPDATE_MESH_HPP
#define CODE_UPDATE_MESH_HPP

#include "boundary_conditions.hpp"
#include "../Rigid_body/body.hpp"

bool fluid_moves(const double t) {
    return true;
}

vec3 global_forces(const double t) {
    return {0.00001*sin(t), 0, 0};
    //return vec3(0);
}

//updating mesh also requires updating the vectors
// - have to extrapolate p, v_n but also v_n1 because some equations require it
template <unsigned N, unsigned M, unsigned P>
void update_mesh(boundary_conditions<N,M,P> &bc, body *b, big_vec<N,M,P, vec3> &v_n, big_vec<N,M,P, vec3> &v_n1, big_vec<N,M,P, double> &p, const double dt, const double t) {
    const auto dx = bc.vel_bc.dx;
    const auto dy = bc.vel_bc.dy;
    const auto dz = bc.vel_bc.dz;


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
                    if (i > 1 && i <N-1 && j > 1 && j< M-1 && k>1 && k<P-1) {  //if off the boundary
                        if (bc.norms.contains(i,j,k) && bc.norms.normal(i,j,k) != vec3(0)) {    //if boundary point outside of mesh
                            points[forces_counter] = bc.points.get_point(i,j,k);    //forces applied at the mesh boundary
                            //taking the force as pointing against the normal, not sure if this is right
                            if (fluid_moves(t)) {
                                forces[forces_counter++] = p(i,j,k)*dx*dy*dz * - bc.norms.normal(i,j,k)  +
                                        global_forces(t); //P=F/A  =>  F=PA
                            } else {
                                forces[forces_counter++] = global_forces(t);
                            }
                        }
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


    b->update_pos(forces, points, dt);

    //extrapolate must be called before update_mesh_boundary because this requires to old values for the normals
    bc.extrapolate(v_n);
    bc.extrapolate(v_n1);
    bc.extrapolate(p);

    //bc.enforce_velocity_BC(v_n1);
    bc.update_mesh_boundary();

    //updating pressure and wall velocity points
    // - non-wall velocity points and normals have already been updated
    //bc.update_velocity_wall_BC();
    //bc.update_pressure_BC(p);




    //can't think of a better way to make sure that the extrapolation does not affect points that need to have BC enforced
    bc.enforce_velocity_BC(v_n);
    //bc.enforce_pressure_BC(p);
}


#endif //CODE_UPDATE_MESH_HPP
