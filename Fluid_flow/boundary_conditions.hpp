//
// Created by jacob on 24/8/21.
//

#ifndef CODE_BOUNDARY_CONDITIONS_HPP
#define CODE_BOUNDARY_CONDITIONS_HPP

#include "../MyMath/boundary.hpp"
#include "../Rigid_body/body.hpp"
#include "../Rigid_body/triangle_mesh.hpp"

template <unsigned N, unsigned M, unsigned P>
struct boundary_conditions {
    unsigned num = 0;   //number of boundary points
    boundary_points<N,M,P> bound;
    boundary_normals<N,M,P> norms;

    big_vec<N,M,P,vec3> vel_bc;
    big_vec<N,M,P,double> p_bc;

    triangle_mesh tm;



    boundary_conditions(const mesh *m, const double dx, const double dy, const double dz) : tm(m) {
        std::cerr << "Remember to change velocity BC back to 0\n";

        set_wall_points();
        norms = boundary_normals<N,M,P>(num);
        create_wall_normals();

        p_bc = big_vec<N,M,P,double>(dx, dy, dz, &bound);
        vel_bc = big_vec<N,M,P,vec3>(dx, dy, dz, &bound);

#ifndef NDEBUG
        DEBUG_check_normal_for_all_boundary_points();
#endif
    }

    void update_pressure_BC();
    void update_velocity_BC();

    void enforce_velocity_BC(big_vec<N,M,P,vec3> &v) {
        for (unsigned i = 0; i <= N; i++) {
            for (unsigned j = 0; j <= M; j++) {
                for (unsigned k = 0; k <= P; k++) {
                    if (bound.is_boundary(i,j,k)) {
                        v.add_elm(i,j,k, vel_bc(i,j,k));
                    }

                }
            }
        }
    }

    void enforce_pressure_BC(big_vec<N,M,P,double> &p) {
        for (unsigned i = 0; i <= N; i++) {
            for (unsigned j = 0; j <= M; j++) {
                for (unsigned k = 0; k <= P; k++) {
                    if (bound.is_boundary(i,j,k)) {
                        p(i,j,k) = p_bc(i,j,k);
                    }

                }
            }
        }
    }

private:
    void set_wall_points();
    void create_wall_normals();
    void DEBUG_check_normal_for_all_boundary_points() noexcept;
};


template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::set_wall_points() {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(i,0,k).has_down = false;
            bound(i,M,k).has_up = false;
            ++num;
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            bound(i,j,0).has_front = false;
            bound(i,j,P).has_back = false;
            ++num;
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(0,j,k).has_left = false;
            bound(N,j,k).has_right = false;
            ++num;
        }
    }
}



template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::update_pressure_BC() {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                if (p_bc.is_boundary(i,j,k)) {

                    const auto norm = norms.normal(i,j,k);

                    const auto nx = norm.x();
                    const auto ny = norm.y();
                    const auto nz = norm.z();

                    const auto dx = p_bc.dx;
                    const auto dy = p_bc.dy;
                    const auto dz = p_bc.dz;

                    //picking the direction
                    unsigned big_dir = 0;
                    if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                        big_dir = 1;
                    }
                    if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                        big_dir = 2;
                    }



                    if (big_dir == 0) { //x direction biggest
                        if (!p_bc.has_left(i,j,k)) {   //forward difference
                            p_bc(i,j,k) = ny/nx* smart_deriv<0,1,0>(p_bc, i,j,k)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p_bc, i,j,k)*2*dx/3 - p_bc(i+2,j,k)/3 + 4*p_bc(i+1,j,k)/3;
                        } else if (!p_bc.has_right(i,j,k)) {   //backward difference
                            p_bc(i,j,k) =-ny/nx* smart_deriv<0,1,0>(p_bc, i,j,k)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p_bc, i,j,k)*2*dx/3 - p_bc(i-2,j,k)/3 + 4*p_bc(i-1,j,k)/3;
                        } else {    //central difference
                            p_bc(i+1,j,k) = -ny/nx * smart_deriv<0,1,0>(p_bc, i,j,k)*dx - nz/nx * smart_deriv<0,0,1>(p_bc, i,j,k)*dx + p_bc(i-1,j,k);
                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p_bc.has_down(i, j, k)) { //forward difference
                            p_bc(i,j,k) = nx/ny* smart_deriv<1,0,0>(p_bc, i,j,k)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p_bc, i,j,k)*2*dy/3 - p_bc(i,j+2,k)/3 + 4*p_bc(i,j+1,k)/3;
                        } else if (!p_bc.has_up(i, j, k)) {  //backward difference
                            p_bc(i,j,k) =-nx/ny* smart_deriv<1,0,0>(p_bc, i,j,k)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p_bc, i,j,k)*2*dy/3 - p_bc(i,j-2,k)/3 + 4*p_bc(i,j-1,k)/3;
                        } else {  //central difference
                            p_bc(i,j+1,k) = -nx/ny * smart_deriv<1,0,0>(p_bc, i,j,k)*dy - nz/ny * smart_deriv<0,0,1>(p_bc, i,j,k)*dy + p_bc(i,j-1,k);
                        }
                    } else {    //z direction
#ifndef NDEBUG
                        if (big_dir > 2) {
                            std::cerr << "the biggest direction cannot be larger than 2\n";
                        }
#endif
                        if (!p_bc.has_front(i,j,k)) {  //forward difference
                            p_bc(i,j,k) = ny/nz* smart_deriv<0,1,0>(p_bc, i,j,k)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p_bc, i,j,k)*2*dz/3 - p_bc(i,j,k+2)/3 + 4*p_bc(i,j,k+1)/3;
                        } else if (!p_bc.has_back(i,j,k)) { //backward difference
                            p_bc(i,j,k) =-ny/nz* smart_deriv<0,1,0>(p_bc, i,j,k)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p_bc, i,j,k)*2*dz/3 - p_bc(i,j,k-2)/3 + 4*p_bc(i,j,k-1)/3;
                        } else {
                            p_bc(i,j,k+1) = -ny/nz * smart_deriv<0,1,0>(p_bc, i,j,k)*dz - nx/nz * smart_deriv<1,0,0>(p_bc, i,j,k)*dz + p_bc(i,j,k-1);
                        }
                    }

                }

            }
        }
    }
}


template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::update_velocity_BC()  {
    //setting the walls
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            vel_bc.add_elm(i,j,0, 0,0,0);
            vel_bc.add_elm(i,j,P, 0,0,0);
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            vel_bc.add_elm(0,j,k, 0,0,0);
            vel_bc.add_elm(N,j,k, 0,0,0);
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            vel_bc.add_elm(i, 0, k, 0.1, 0, 0);
            //vel_bc.add_elm(i, 0, k, 0, 0, 0);
            vel_bc.add_elm(i,M,k, 0,0,0);
        }
    }

}

template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::create_wall_normals() {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            norms.add_point(i,0,k, vec3(0,1,0));
            norms.add_point(i,M,k, vec3(0,-1,0));
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            norms.add_point(i,j,0, vec3(0,0,1));
            norms.add_point(i,j,P, vec3(0,0,-1));
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            norms.add_point(0,j,k,vec3(1,0,0));
            norms.add_point(N,j,k,vec3(-1,0,0));
        }
    }
}

template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::DEBUG_check_normal_for_all_boundary_points() noexcept {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <=M; j++) {
            for (unsigned k = 0; k<=P; k++) {
                if (bound.is_boundary(i,j,k)) {
                    if (!norms.contains(i,j,k)) {
                        std::cerr << "Boundary at i=" << i << " j=" << j << " k=" << k << "does not have a normal\n";
                    }
                }
            }
        }
    }
}

#endif //CODE_BOUNDARY_CONDITIONS_HPP
