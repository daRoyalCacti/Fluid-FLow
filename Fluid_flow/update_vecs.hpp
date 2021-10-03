//
// Created by jacob on 15/9/21.
//

#ifndef CODE_UPDATE_VECS_HPP
#define CODE_UPDATE_VECS_HPP

#include "boundary_conditions.hpp"
#include "../MyMath/big_vec.hpp"


//v is vector to enforce the conditions on
void enforce_velocity_BC(const boundary_conditions &BC,  big_vec_v &v) {
    //std::cerr << "testing\n";
    for (unsigned i = 0; i < v.size(); i++) {
        //std::cerr << "\t" << i << "/" << v.size() << "\t" << v.is_boundary(i) << "\n";
        if (v.is_boundary(i)) {
            //if (BC.v_points.get_vel(i) == vec3(0.75,0,0)) {
            //    std::cerr << "whasdf\n";
            //}
            v.add_elm(i, BC.v_points.get_vel(i));
#ifndef NDEBUG
            if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())) {
                std::cerr << "velocity boundary condition returned an infinite value\n";
            }
#endif
        }

    }


}



void update_pressure_BC(const boundary_conditions &BC, big_vec_d &p) {
        for (unsigned i = 0; i < p.size(); i++) {
            if (p.is_boundary(i)) {
#ifndef NDEBUG
                if (!BC.norms.contains(i)) {
                    std::cerr << "norms does not contain " << i << ")\n";
                }
#endif

                const auto norm = BC.norms.normal(i);

#ifndef NDEBUG
                if (norm == vec3(0) ) {
                    std::cerr << "normal vector is the 0 vector. This should never happen!\n";
                }
#endif

                const auto nx = norm.x();
                const auto ny = norm.y();
                const auto nz = norm.z();

                const auto dx = p.dx(i);
                const auto dy = p.dy(i);
                const auto dz = p.dz(i);

                //picking the direction
                unsigned big_dir = 0;
                if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                    big_dir = 1;
                }
                if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                    big_dir = 2;
                }



                if (big_dir == 0) { //x direction biggest
                    if (!p.has_right(i)) {   //backward difference
                        p(i) =-ny/nx* smart_deriv<0,1,0>(p, i)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i)*2*dx/3 - p.move(i,-2,0,0)/3 + 4*p.move(i,-1,0,0)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else {   //forward difference
                        p(i) = ny/nx* smart_deriv<0,1,0>(p, i)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i)*2*dx/3 - p.move(i,2,0,0)/3 + 4*p.move(i,1,0,0)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    }
                } else if (big_dir == 1) {  //y direction biggest
                    if (!p.has_up(i)) {  //backward difference
                        p(i) =-nx/ny* smart_deriv<1,0,0>(p, i)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i)*2*dy/3 - p.move(i,0,-2,0)/3 + 4*p.move(i,0,-1,0)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else { //forward difference
                        p(i) = nx/ny* smart_deriv<1,0,0>(p, i)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i)*2*dy/3 - p.move(i, 0,2,0)/3 + 4*p.move(i,0,1,0)/3;
 #ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    }
                } else {    //z direction
#ifndef NDEBUG
                    if (big_dir > 2) {
                        std::cerr << "the biggest direction cannot be larger than 2\n";
                    }
#endif
                    if (!p.has_back(i)) { //backward difference
                        p(i) =-ny/nz* smart_deriv<0,1,0>(p, i)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i)*2*dz/3 - p.move(i,0,0,-2)/3 + 4*p.move(i,0,0,-1)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else {  //forward difference
                        p(i) = ny/nz* smart_deriv<0,1,0>(p, i)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i)*2*dz/3 - p.move(i,0,0,2)/3 + 4*p.move(i,0,0,1)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    }


                }

            }

        }



}

#endif //CODE_UPDATE_VECS_HPP
