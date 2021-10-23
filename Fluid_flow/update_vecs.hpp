//
// Created by jacob on 15/9/21.
//

#ifndef CODE_UPDATE_VECS_HPP
#define CODE_UPDATE_VECS_HPP

#ifndef NDEBUG
#define UPDATE_VECS_CHECK_RESULTS_LOG
#endif

#include "boundary_conditions.hpp"
#include "../MyMath/big_vec.hpp"

template <bool err>
bool enforce_velocity_correction_BC(const boundary_conditions &BC,  big_vec_v &v, const double accuracy_percent = 1.0) {
    for (unsigned i = 0; i < v.size(); i++) {
        if (v.is_boundary(i) ) {
            if (!v.g->off_walls(i)) {
                const auto old = v(i);

#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                if (!BC.norms.contains(i)) {
                    std::cerr << "norms does not contain " << i << ")\n";
                    return false;
                }
#endif

                const auto norm = BC.norms.normal(i);

#ifndef NDEBUG
                if (norm == vec3(0) ) {
                    std::cerr << "normal vector is the 0 vector. This should never happen!\n";
                    return false;
                }
#endif

                v.add_elm(i, 4/3.0*v.move(i, norm) - 1/3.0*v.move(i,2*norm) );

                const bool is_good_x = old.x() > 1e-3 && std::abs( (old.x() - v(i).x()) /v(i).x() *100) > accuracy_percent;
                const bool is_good_y = old.y() > 1e-3 && std::abs( (old.y() - v(i).y()) /v(i).y() *100) > accuracy_percent;
                const bool is_good_z = old.z() > 1e-3 && std::abs( (old.z() - v(i).z()) /v(i).z() *100) > accuracy_percent;
                const bool is_good = is_good_x && is_good_y && is_good_z;

#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())  ) {
                    std::cerr << "velocity boundary condition returned an infinite value\n";
                }
                if constexpr (err) {
                    if (is_good) {
                        std::cerr << "Too much correction on the velocity correction at the walls\n";
                        std::cerr << "\tcorrection : " << ((old-v(i))/v(i)).length()*100 << "%\n";
                        std::cerr << "\tat index " << i << "\n";
                        std::cerr << "\told velocity = (" << old << ")\t new velocity = (" << v(i) << ")\n";
                    }
                }
#endif
                if constexpr (err) {
                    if (is_good) {
                        return false;
                    }
                }
                if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())  ) {
                    return false;
                }

            }
        }


    }

    return true;
}


//v is vector to enforce the conditions on
template <bool err>
bool enforce_velocity_BC(const boundary_conditions &BC,  big_vec_v &v, const double accuracy_percent = 1.0) {
    for (unsigned i = 0; i < v.size(); i++) {
        if (v.is_boundary(i) ) {
            const auto old = v(i);
            if (v.g->off_walls(i)) {
                v.add_elm(i, BC.v_points.get_vel(i));
#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())) {
                    std::cerr << "velocity boundary condition returned an infinite value\n";
                    return false;
                }
                if constexpr (err) {
                    if ((old-v(i)).length()/v(i).length()*100 > accuracy_percent) {
                        std::cerr << "Too much correction on the velocity at the mesh\n";
                        std::cerr << "\tcorrection : " << (old-v(i)).length()/v(i).length()*100 << "%\n";
                        std::cerr << "\tat index " << i << "\n";
                        std::cerr << "\told velocity = (" << old << ")\t new velocity = (" << v(i) << ")\n";
                    }
                }

#endif

            } else {

#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                if (!BC.norms.contains(i)) {
                    std::cerr << "norms does not contain " << i << ")\n";
                    return false;
                }
#endif
                const auto norm = BC.norms.normal(i);
                const auto old_deriv = norm.x()*smart_deriv_old<1,0,0>(v, i) + norm.y()*smart_deriv_old<0,1,0>(v, i) + norm.z()*smart_deriv_old<0,0,1>(v, i);


#ifndef NDEBUG
                if (norm == vec3(0) ) {
                    std::cerr << "normal vector is the 0 vector. This should never happen!\n";
                    return false;
                }
#endif

                v.add_elm(i, 4/3.0*v.move(i, norm) - 1/3.0*v.move(i,2*norm) );



                const bool is_good_x = old.x() > 1e-3 && std::abs( (old.x() - v(i).x()) /v(i).x() *100) > accuracy_percent;
                const bool is_good_y = old.y() > 1e-3 && std::abs( (old.y() - v(i).y()) /v(i).y() *100) > accuracy_percent;
                const bool is_good_z = old.z() > 1e-3 && std::abs( (old.z() - v(i).z()) /v(i).z() *100) > accuracy_percent;
                const bool is_good = is_good_x && is_good_y && is_good_z;

#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())  ) {
                    std::cerr << "velocity boundary condition returned an infinite value\n";
                }
                if constexpr (err) {
                    if (is_good) {
                        std::cerr << "Too much correction on the velocity at the walls\n";

                        std::cerr << "\tThis was done using : ";
                        if (!v.has_right(i)) {
                            std::cerr << "backwards in x,  ";
                        } else {
                            std::cerr << "forwards in x,  ";
                        }
                        if (!v.has_up(i)) {
                            std::cerr << "backwards in y,  ";
                        } else {
                            std::cerr << "forwards in y,  ";
                        }
                        if (!v.has_back(i)) {
                            std::cerr << "backwards in z\n";
                        } else {
                            std::cerr << "forwards in z\n";
                        }

                        std::cerr << "\tcorrection : " << ((old-v(i))/v(i)).length()*100 << "%\n";
                        std::cerr << "\tat index " << i << "\n";
                        std::cerr << "\told velocity = (" << old << ")\t new velocity = (" << v(i) << ")\n";
                        std::cerr << "\tnormal vector = (" << norm << ")\n";


                        std::cerr << "\told derivative (" << old_deriv << ")\tnew derivative (" << norm.x()*smart_deriv_old<1,0,0>(v, i) + norm.y()*smart_deriv_old<0,1,0>(v, i) + norm.z()*smart_deriv_old<0,0,1>(v, i) << ")\n";
                        std::cerr << "\tThis is at (" << BC.global_grid.get_plot_pos(i) << ")\n";
                        std::cerr << "\tl:" << v.has_left(i) << " r:" << v.has_right(i) << " u:" << v.has_up(i) << " d:" << v.has_down(i) << " f:" << v.has_front(i) << " b:" << v.has_back(i) << "\n";

                    }
                }
#endif

            }

            const bool is_good_x = old.x() > 1e-3 && std::abs( (old.x() - v(i).x()) /v(i).x() *100) > accuracy_percent;
            const bool is_good_y = old.y() > 1e-3 && std::abs( (old.y() - v(i).y()) /v(i).y() *100) > accuracy_percent;
            const bool is_good_z = old.z() > 1e-3 && std::abs( (old.z() - v(i).z()) /v(i).z() *100) > accuracy_percent;
            const bool is_good = is_good_x && is_good_y && is_good_z;

            if constexpr (err) {
                if (is_good) {
                    return false;
                }
            }
            if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())  ) {
                return false;
            }
        }


    }

    return true;
}



template <bool err>
bool update_pressure_BC(const boundary_conditions &BC, big_vec_d &p, const double accuracy_percent = 1.0) {
        for (unsigned i = 0; i < p.size(); i++) {
            if (p.is_boundary(i)) {
                const auto old = p(i);

                if (!BC.global_grid.off_walls(i)) {
                    p(i) = 0;
#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                    if constexpr (err) {
                        if (old > 1e-4 && std::abs(old-p(i))/std::abs(p(i))*100 > accuracy_percent) {
                            std::cerr << "solution to pressure is not accurate\n";

                            std::cerr << "\tcorrection : " << std::abs(old-p(i))/std::abs(p(i))*100 << "%\n";
                            std::cerr << "\tat index " << i << "\n";
                            std::cerr << "\told pressure = " << old << "\t new pressure = " << p(i) << "\n";
                            std::cerr << "\tbad at edge point\n";
                        }
                    }
#endif
                } else {

#ifndef NDEBUG
                    if (!BC.norms.contains(i)) {
                        std::cerr << "norms does not contain " << i << ")\n";
                        return false;
                    }
#endif

                    const auto norm = BC.norms.normal(i);

#ifndef NDEBUG
                    if (norm == vec3(0) ) {
                        std::cerr << "normal vector is the 0 vector. This should never happen!\n";
                        return false;
                    }
#endif

                    const auto nx = norm.x();
                    const auto ny = norm.y();
                    const auto nz = norm.z();

                    const auto dx = p.dx(i);
                    const auto dy = p.dy(i);
                    const auto dz = p.dz(i);

                    const bool cx = p.can_move(i, -1,0,0) && p.can_move(i, 1,0,0);
                    const bool cy = p.can_move(i, 0,-1,0) && p.can_move(i, 0,1,0);
                    const bool cz = p.can_move(i, 0,0,-1) && p.can_move(i, 0,0,1);

                    const bool fx = p.can_move(i, 2,0,0);
                    const bool fy = p.can_move(i, 0,2,0);
                    const bool fz = p.can_move(i, 0,0,2);

                    const bool bx = p.can_move(i, -2,0,0);
                    const bool by = p.can_move(i, 0,-2,0);
                    const bool bz = p.can_move(i, 0,0,-2);

                    std::string_view err_code;

                    //cc not possible because have boundary point
                    if (cx & cy & fz) { //ccf
                        err_code = "ccf";
                        p(i) = 2*dz/3 * ( nx/nz * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)   + ny/nz * (p.move(i,0,1,0) - p.move(i,0,-1,0)) / (2*dy)  ) + 4/3.0*p.move(i,0,0,1) - 1/3.0*p.move(i,0,0,2);
                    } else if (cx & fy & cz) {  //cfc
                        err_code = "cfc";
                        p(i) = 2*dy/3 * ( nx/ny * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)   + nz/ny * (p.move(i,0,0,1) - p.move(i,0,0,-1)) / (2*dz)  ) + 4/3.0*p.move(i,0,1,0) - 1/3.0*p.move(i,0,2,0);
                    } else if (fx & cy & cz) {  //fcc
                        err_code = "fcc";
                        p(i) = 2*dx/3 * ( ny/nx * (p.move(i,0,1,0) - p.move(i,0,-1,0))/(2*dy)   + nz/nx * (p.move(i,0,0,1) - p.move(i,0,0,-1)) / (2*dz)  ) + 4/3.0*p.move(i,1,0,0) - 1/3.0*p.move(i,2,0,0);

                    } else if (cx & cy & bz) {  //ccb
                        err_code = "ccb";
                        p(i) = 2*dz/3 * ( -nx/nz * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)   - ny/nz * (p.move(i,0,1,0) - p.move(i,0,-1,0)) / (2*dy)  ) + 4/3.0*p.move(i,0,0,-1) - 1/3.0*p.move(i,0,0,-2);
                    } else if (cx & by & cz) {  //cbc
                        err_code = "cbc";
                        p(i) = 2*dy/3 * ( -nx/ny * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)   - nz/ny * (p.move(i,0,0,1) - p.move(i,0,0,-1)) / (2*dz)  ) + 4/3.0*p.move(i,0,-1,0) - 1/3.0*p.move(i,0,-2,0);
                    } else if (bx & cy & cz) {  //bcc
                        err_code = "bcc";
                        p(i) = 2*dx/3 * ( -ny/nx * (p.move(i,0,1,0) - p.move(i,0,-1,0))/(2*dy)   - nz/nx * (p.move(i,0,0,1) - p.move(i,0,0,-1)) / (2*dz)  ) + 4/3.0*p.move(i,-1,0,0) - 1/3.0*p.move(i,-2,0,0);

                    } else if (cx & fy & fz) {  //cff
                        err_code = "cff";
                        p(i) = 1 / ( 3*nz/(2*dz)  + 3*ny/(2*dy) ) * (  nx * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)  +  ny * (-p.move(i,0,2,0) + 4*p.move(i,0,1,0))/(2*dy)  + nz * (-p.move(i,0,0,2)+4*p.move(i,0,0,1))/(2*dz)  );
                    } else if (fx & cy & fz) {  //fcf
                        err_code = "fcf";
                        p(i) = 1 / ( 3*nx/(2*dx)  + 3*nz/(2*dz) ) * (  ny * (p.move(i,0,1,0) - p.move(i,0,-1,0))/(2*dy)  +  nx * (-p.move(i,2,0,0) + 4*p.move(i,1,0,0))/(2*dx)  + nz * (-p.move(i,0,0,2)+4*p.move(i,0,0,1))/(2*dz)  );
                    } else if (fx & fy & cz) {  //ffc
                        err_code = "ffc";
                        p(i) = 1 / ( 3*nx/(2*dx)  + 3*ny/(2*dy) ) * (  nz * (p.move(i,0,0,1) - p.move(i,0,0,-1))/(2*dz)  +  ny * (-p.move(i,0,2,0) + 4*p.move(i,0,1,0))/(2*dy)  + nx * (-p.move(i,2,0,0)+4*p.move(i,1,0,0))/(2*dx)  );

                    } else if (cx & fy & bz) {  //cfb
                        err_code = "cfb";
                        p(i) = 1 / ( 3*nz/(2*dz)  - 3*ny/(2*dy) )  *  (  -nx * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)  -  ny * (-p.move(i,0,2,0)+4*p.move(i,0,1,0))/(2*dy) - nz*(-4*p.move(i,0,0,-1)+p.move(i,0,0,-2))/(2*dz)  );
                    } else if (fx & cy & bz) {  //fcb
                        err_code = "fcb";
                        p(i) = 1 / ( 3*nz/(2*dz)  - 3*nx/(2*dx) )  *  (  -ny * (p.move(i,0,1,0) - p.move(i,0,-1,0))/(2*dy)  -  nx * (-p.move(i,2,0,0)+4*p.move(i,1,0,0))/(2*dx) - nz*(-4*p.move(i,0,0,-1)+p.move(i,0,0,-2))/(2*dz)  );
                    } else if (fx & by & cz) {  //fbc
                        err_code = "fbc";
                        p(i) = 1 / ( 3*ny/(2*dy)  - 3*nx/(2*dx) )  *  (  -nz * (p.move(i,0,0,1) - p.move(i,0,0,-1))/(2*dz)  -  nx * (-p.move(i,2,0,0)+4*p.move(i,1,0,0))/(2*dx) - ny*(-4*p.move(i,0,-1,0)+p.move(i,0,-2,0))/(2*dy)  );
                    } else if (cx & by & fz) {  //cbf
                        err_code = "cbf";
                        p(i) = 1 / ( 3*ny/(2*dy)  - 3*nz/(2*dz) )  *  (  -nx * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx)  -  nz * (-p.move(i,0,0,2)+4*p.move(i,0,0,1))/(2*dz) - ny*(-4*p.move(i,0,-1,0)+p.move(i,0,-2,0))/(2*dy)  );
                    } else if (bx & cy & fz) {  //bcf
                        err_code = "bcf";
                        p(i) = 1 / ( 3*nx/(2*dx)  - 3*nz/(2*dz) )  *  (  -ny * (p.move(i,0,1,0) - p.move(i,0,-1,0))/(2*dy)  -  nz * (-p.move(i,0,0,2)+4*p.move(i,0,0,1))/(2*dz) - nx*(-4*p.move(i,-1,0,0)+p.move(i,-2,0,0))/(2*dx)  );
                    } else if (bx & fy & cz) {  //bfc
                        err_code = "bfc";
                        p(i) = 1 / ( 3*nx/(2*dx)  - 3*ny/(2*dy) )  *  (  -nz * (p.move(i,0,0,1) - p.move(i,0,0,-1))/(2*dz)  -  ny * (-p.move(i,0,2,0)+4*p.move(i,0,1,0))/(2*dy) - nx*(-4*p.move(i,-1,0,0)+p.move(i,-2,0,0))/(2*dx)  );

                    } else if (cx & by & bz) {  //cbb
                        err_code = "cbb";
                        p(i) = 1 / ( 3*ny/(2*dy) + 3*nz/(2*nz) )  *  (  -nx * (p.move(i,1,0,0) - p.move(i,-1,0,0))/(2*dx) - ny * (-4*p.move(i,0,-1,0) + p.move(i,0,0,-2))/(2*dy) - nz*(-4*p.move(i,0,0,-1)+p.move(i,0,0,-2))/(2*dz) );
                    } else if (bx & cy & bz) {  //bcb
                        err_code = "bcb";
                        p(i) = 1 / ( 3*nx/(2*dx) + 3*nz/(2*nz) )*  (  -ny * (p.move(i,0,1,0) - p.move(i,0,-1,0))/(2*dy) - nx * (-4*p.move(i,-1,0,0) + p.move(i,0,-2,0))/(2*dx) - nz*(-4*p.move(i,0,0,-1)+p.move(i,0,0,-2))/(2*dz) );
                    } else if (bx & by & cz) {  //bbc
                        err_code = "bbc";
                        p(i) = 1 / ( 3*ny/(2*dy) + 3*nx/(2*nx) )*  (  -nz * (p.move(i,0,0,1) - p.move(i,0,0,-1))/(2*dz) - ny * (-4*p.move(i,0,-1,0) + p.move(i,0,0,-2))/(2*dy) - nx*(-4*p.move(i,-1,0,0)+p.move(i,-2,0,0))/(2*dx) );

                    } else if (fx & fy & fz) {  //fff
                        err_code = "fff";
                        p(i) = 1/ ( 3*nx/(2*dx) + 3*ny/(2*ny) + 3*nz/(2*nz) ) * (  nx*(-p.move(i,2,0,0)+4*p.move(i,1,0,0))/(2*dx)  +  ny*(-p.move(i,0,2,0)+4*p.move(i,0,1,0))/(2*dy)  +  nz*(-p.move(i,0,0,2)+4*p.move(i,0,0,1))/(2*dz)  );
                    } else if (bx & by & bz) {  //bbb
                        err_code = "bbb";
                        p(i) = 1/ ( 3*nx/(2*dx) + 3*ny/(2*ny) + 3*nz/(2*nz) ) * (  nx*(-p.move(i,-2,0,0)+4*p.move(i,-1,0,0))/(2*dx)  +  ny*(-p.move(i,0,-2,0)+4*p.move(i,0,-1,0))/(2*dy)  +  nz*(-p.move(i,0,0,-2)+4*p.move(i,0,0,-1))/(2*dz)  );

                    } else if (fx & fy & bz) {  //ffb
                        err_code = "ffb";
                        p(i) = 1/( -3*nx/(2*dx) - 3*ny/(2*dy) + 3*nz/(2*dz) )  *  (   -nx*(-p.move(i, 2, 0,0) + 4*p.move(i,1,0,0) )/(2*dx)  -  ny * (-p.move(i, 0, 2, 0) + 4*p.move(i, 0, 1,0))/(2*dy)  +  nz * (-p.move(i, 0, 0, -2) + 4*p.move(i, 0, 0,-1))/(2*dz) );
                    } else if (fx & by & fz) {  //fbf
                        err_code = "fbf";
                        p(i) = 1/( -3*nx/(2*dx) + 3*ny/(2*dy) - 3*nz/(2*dz) ) *  (   -nx*(-p.move(i, 2, 0,0) + 4*p.move(i,1,0,0) )/(2*dx)  +  ny * (-p.move(i, 0, -2, 0) - 4*p.move(i, 0, -1,0))/(2*dy)  -  nz * (-p.move(i, 0, 0, 2) + 4*p.move(i, 0, 0,1))/(2*dz) );
                    } else if (bx & fy & fz) {  //bff
                        err_code = "bff";
                        p(i) = 1/( 3*nx/(2*dx) - 3*ny/(2*dy) - 3*nz/(2*dz) ) *  (   nx*(-p.move(i, -2, 0,0) + 4*p.move(i,-1,0,0) )/(2*dx)  -  ny * (-p.move(i, 0, 2, 0) - 4*p.move(i, 0, 1,0))/(2*dy)  +  nz * (-p.move(i, 0, 0, 2) + 4*p.move(i, 0, 0,1))/(2*dz) );

                    } else if (bx & by & fz) {  //bbf
                        err_code = "bbf";
                        p(i) = 1/( 3*nx/(2*dx) + 3*ny/(2*dy) - 3*nz/(2*dz) )  *  (   nx*(-p.move(i, -2, 0,0) + 4*p.move(i,-1,0,0) )/(2*dx)  +  ny * (-p.move(i, 0, -2, 0) + 4*p.move(i, 0, -1,0))/(2*dy)  -  nz * (-p.move(i, 0, 0, 2) + 4*p.move(i, 0, 0,1))/(2*dz) );
                    } else if (bx & fy & bz) {  //bfb
                        err_code = "bfb";
                        p(i) = 1/( 3*nx/(2*dx) - 3*ny/(2*dy) + 3*nz/(2*dz) )  *  (   nx*(-p.move(i, -2, 0,0) + 4*p.move(i,-1,0,0) )/(2*dx)  -  ny * (-p.move(i, 0, 2, 0) + 4*p.move(i, 0, 1,0))/(2*dy)  +  nz * (-p.move(i, 0, 0, -2) + 4*p.move(i, 0, 0,-1))/(2*dz) );
                    } else if (fx & by & bz) {  //fbb
                        err_code = "fbb";
                        p(i) = 1/( -3*nx/(2*dx) + 3*ny/(2*dy) + 3*nz/(2*dz) )  *  (   -nx*(-p.move(i, 2, 0,0) + 4*p.move(i,1,0,0) )/(2*dx)  +  ny * (-p.move(i, 0, -2, 0) + 4*p.move(i, 0, -1,0))/(2*dy)  +  nz * (-p.move(i, 0, 0, -2) + 4*p.move(i, 0, 0,-1))/(2*dz) );
                    } else {
                        std::cerr << "enforcing pressure boundary conditions is not possible.\n";
                    }
#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                    if (!std::isfinite(p(i))) {
                        std::cerr << "pressure boundary condition returned an infinite value\n";
                    }
                    if constexpr (err) {
                        if (old > 1e-4 && std::abs(old-p(i))/std::abs(p(i))*100 > accuracy_percent) {
                            std::cerr << "solution to pressure is not accurate\n";

                            std::cerr << "\tcorrection : " << std::abs(old-p(i))/std::abs(p(i))*100 << "%\n";
                            std::cerr << "\tat index " << i << "\n";
                            std::cerr << "\told pressure = " << old << "\t new pressure = " << p(i) << "\n";
                            std::cerr << "\tnormal vector = (" << norm << ")\n";

                            std::cerr << "\tThis was done using : " << err_code << "\n";
                        }
                    }
#endif
                }

                if constexpr (err) {
                    if (old > 1e-4 && std::abs(old-p(i))/std::abs(p(i))*100 > accuracy_percent) {
                        return false;
                    }
                }

                if (!std::isfinite(p(i))) {
                    return false;
                }


            }

        }


    return true;
}

#endif //CODE_UPDATE_VECS_HPP
