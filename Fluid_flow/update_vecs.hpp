//
// Created by jacob on 15/9/21.
//

#ifndef CODE_UPDATE_VECS_HPP
#define CODE_UPDATE_VECS_HPP

#include "boundary_conditions.hpp"
#include "../MyMath/big_vec.hpp"

void update_pressure_BC(const boundary_conditions &BC, big_vec_d &p) {
        for (unsigned i = 0; i < p.size(); i++) {
            if (p.is_boundary(i)) {
#ifndef NDEBUG
                if (!BC.norms.contains(i)) {
                    std::cerr << "norms does not contain " << i << ")\n";
                }
#endif

const auto norm = BC.norms.normal(i);

                //this occurs when inside a boundary
                if (norm == vec3(0) ) {
                    p(i) = 0;
                    continue;   //rest of the code will only error
                }

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
                    if (!p.has_left(i)) {   //forward difference
                        p(i) = ny/nx* smart_deriv<0,1,0>(p, i)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i)*2*dx/3 - p.move(i,2,0,0)/3 + 4*p.move(i,1,0,0)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else if (!p.has_right(i)) {   //backward difference
                        p(i) =-ny/nx* smart_deriv<0,1,0>(p, i)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i)*2*dx/3 - p.move(i,-2,0,0)/3 + 4*p.move(i,-1,0,0)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else {    //central difference
                        p( p.get_move_ind(i,1,0,0) ) = -ny/nx * smart_deriv<0,1,0>(p, i)*dx - nz/nx * smart_deriv<0,0,1>(p, i)*dx + p.move(i,-1,0,0);
#ifndef NDEBUG
                        if (!std::isfinite(p.move(i,1,0,0))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    }
                } else if (big_dir == 1) {  //y direction biggest
                    if (!p.has_down(i)) { //forward difference
                        p(i) = nx/ny* smart_deriv<1,0,0>(p, i)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i)*2*dy/3 - p.move(i, 0,2,0)/3 + 4*p.move(i,0,1,0)/3;
 #ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else if (!p.has_up(i)) {  //backward difference
                        p(i) =-nx/ny* smart_deriv<1,0,0>(p, i)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i)*2*dy/3 - p.move(i,0,-2,0)/3 + 4*p.move(i,0,-1,0)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else {  //central difference
                        p( p.get_move_ind(i, 0,1,0) ) = -nx/ny * smart_deriv<1,0,0>(p, i)*dy - nz/ny * smart_deriv<0,0,1>(p, i)*dy + p.move(i,0,-1,0);
#ifndef NDEBUG
                        if (!std::isfinite(p.move(i,0,1,0))) {
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
if (!p.has_front(i)) {  //forward difference
                        p(i) = ny/nz* smart_deriv<0,1,0>(p, i)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i)*2*dz/3 - p.move(i,0,0,2)/3 + 4*p.move(i,0,0,1)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else if (!p.has_back(i)) { //backward difference
                        p(i) =-ny/nz* smart_deriv<0,1,0>(p, i)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i)*2*dz/3 - p.move(i,0,0,-2)/3 + 4*p.move(i,0,0,-1)/3;
#ifndef NDEBUG
                        if (!std::isfinite(p(i))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    } else {
                        p(p.get_move_ind(i,0,0,1)) = -ny/nz * smart_deriv<0,1,0>(p, i)*dz - nx/nz * smart_deriv<1,0,0>(p, i)*dz + p.move(i,0,0,-1);
#ifndef NDEBUG
                        if (!std::isfinite(p.move(i,0,0,1))) {
                            std::cerr << "pressure boundary condition returned an infinite value\n";
                        }
#endif
                    }


                }

            }

        }



}

#endif //CODE_UPDATE_VECS_HPP
