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
    //std::cerr << "testing\n";
    for (unsigned i = 0; i < v.size(); i++) {
        //std::cerr << "\t" << i << "/" << v.size() << "\t" << v.is_boundary(i) << "\n";
        if (v.is_boundary(i) ) {
            if (!v.g->off_walls(i)) {
                const auto old = v(i);

#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
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

                const auto dx = v.dx(i);
                const auto dy = v.dy(i);
                const auto dz = v.dz(i);

                //picking the direction
                unsigned big_dir = 0;
                if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                    big_dir = 1;
                }
                if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                    big_dir = 2;
                }

                if (big_dir == 0) { //x direction biggest
                    if (!v.has_right(i)) {   //backward difference
                        v.add_elm(i, -ny/nx* smart_deriv<0,1,0>(v, i)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(v, i)*2*dx/3 - v.move(i,-2,0,0)/3 + 4*v.move(i,-1,0,0)/3);
                    } else {   //forward difference
                        v.add_elm(i, ny/nx* smart_deriv<0,1,0>(v, i)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(v, i)*2*dx/3 - v.move(i,2,0,0)/3 + 4*v.move(i,1,0,0)/3);
                    }
                } else if (big_dir == 1) {  //y direction biggest
                    if (!v.has_up(i)) {  //backward difference
                        v.add_elm(i, -nx/ny* smart_deriv<1,0,0>(v, i)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(v, i)*2*dy/3 - v.move(i,0,-2,0)/3 + 4*v.move(i,0,-1,0)/3);
                    } else { //forward difference
                        v.add_elm(i, nx/ny* smart_deriv<1,0,0>(v, i)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(v, i)*2*dy/3 - v.move(i, 0,2,0)/3 + 4*v.move(i,0,1,0)/3);
                    }
                } else {    //z direction
#ifndef NDEBUG
                    if (big_dir > 2) {
                        std::cerr << "the biggest direction cannot be larger than 2\n";
                    }
#endif
                    if (!v.has_back(i)) { //backward difference
                        v.add_elm(i, -ny/nz* smart_deriv<0,1,0>(v, i)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(v, i)*2*dz/3 - v.move(i,0,0,-2)/3 + 4*v.move(i,0,0,-1)/3);
                    } else {  //forward difference
                        v.add_elm(i, ny/nz* smart_deriv<0,1,0>(v, i)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(v, i)*2*dz/3 - v.move(i,0,0,2)/3 + 4*v.move(i,0,0,1)/3);
                    }

                }
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
                        std::cerr << "\tcorrection : " << (old-v(i)).length()/v(i).length()*100 << "%\n";
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

            }
        }


    }

    return true;
}


//v is vector to enforce the conditions on
template <bool err>
bool enforce_velocity_BC(const boundary_conditions &BC,  big_vec_v &v, const double accuracy_percent = 1.0) {
    //std::cerr << "testing\n";
    for (unsigned i = 0; i < v.size(); i++) {
        //std::cerr << "\t" << i << "/" << v.size() << "\t" << v.is_boundary(i) << "\n";
        if (v.is_boundary(i) ) {
            const auto old = v(i);
            if (v.g->off_walls(i)) {
                v.add_elm(i, BC.v_points.get_vel(i));
#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                if (!std::isfinite(v(i).x()) || !std::isfinite(v(i).y()) || !std::isfinite(v(i).z())) {
                    std::cerr << "velocity boundary condition returned an infinite value\n";
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

                const auto dx = v.dx(i);
                const auto dy = v.dy(i);
                const auto dz = v.dz(i);

                //picking the direction
                unsigned big_dir = 0;
                if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                    big_dir = 1;
                }
                if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                    big_dir = 2;
                }

                if (big_dir == 0) { //x direction biggest
                    if (!v.has_right(i)) {   //backward difference
                        v.add_elm(i, -ny/nx* smart_deriv<0,1,0>(v, i)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(v, i)*2*dx/3 - v.move(i,-2,0,0)/3 + 4*v.move(i,-1,0,0)/3);
                    } else {   //forward difference
                        v.add_elm(i, ny/nx* smart_deriv<0,1,0>(v, i)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(v, i)*2*dx/3 - v.move(i,2,0,0)/3 + 4*v.move(i,1,0,0)/3);
                    }
                } else if (big_dir == 1) {  //y direction biggest
                    if (!v.has_up(i)) {  //backward difference
                        v.add_elm(i, -nx/ny* smart_deriv<1,0,0>(v, i)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(v, i)*2*dy/3 - v.move(i,0,-2,0)/3 + 4*v.move(i,0,-1,0)/3);
                    } else { //forward difference
                        v.add_elm(i, nx/ny* smart_deriv<1,0,0>(v, i)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(v, i)*2*dy/3 - v.move(i, 0,2,0)/3 + 4*v.move(i,0,1,0)/3);
                    }
                } else {    //z direction
#ifndef NDEBUG
                    if (big_dir > 2) {
                        std::cerr << "the biggest direction cannot be larger than 2\n";
                    }
#endif
                    if (!v.has_back(i)) { //backward difference
                        v.add_elm(i, -ny/nz* smart_deriv<0,1,0>(v, i)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(v, i)*2*dz/3 - v.move(i,0,0,-2)/3 + 4*v.move(i,0,0,-1)/3);
                    } else {  //forward difference
                        v.add_elm(i, ny/nz* smart_deriv<0,1,0>(v, i)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(v, i)*2*dz/3 - v.move(i,0,0,2)/3 + 4*v.move(i,0,0,1)/3);
                    }

                }
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
                        if (big_dir == 0) {
                            std::cerr << "\tThis happened for x the largest direction\n";
                        } else if (big_dir == 1) {
                            std::cerr << "\tThis happened for y the largest direction\n";
                        } else {
                            std::cerr << "\tThis happened for z the largest direction\n";
                        }
                        std::cerr << "\tcorrection : " << (old-v(i)).length()/v(i).length()*100 << "%\n";
                        std::cerr << "\tat index " << i << "\n";
                        std::cerr << "\told velocity = (" << old << ")\t new velocity = (" << v(i) << ")\n";
                        std::cerr << "\tnormal vector = (" << norm << ")\n";

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
                            std::cerr << "backwards in z";
                        } else {
                            std::cerr << "forwards in z";
                        }
                        std::cerr << "\n";
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
        }


    }

    return true;
}



/*void update_wall_velocity(boundary_conditions &bc, const big_vec_v &v) {
    const auto dims = bc.global_grid.no_points_unif;
    const auto N = static_cast<unsigned>(dims.x()-1);
    const auto M = static_cast<unsigned>(dims.y()-1);
    const auto P = static_cast<unsigned>(dims.z()-1);

    const auto dx = v.g->dx;
    const auto dy = v.g->dy;
    const auto dz = v.g->dz;

    const auto &vx = v.xv;
    const auto &vy = v.yv;
    const auto &vz = v.zv;


    //setting the walls
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            //need to change v.z
            const auto ind1 = bc.old_new[bc.global_grid.convert_indices_unif(vec3(i,j,0) )];
            //const auto new_vec1 = vec3( dx*smart_deriv<0,1,0>(vy, ind1) + dx*smart_deriv<0,0,1>(vz, ind1) + vz.move(ind1, 1,0,0), vy(ind1), vz(ind1));
            const auto new_vec1 = vec3( vx(ind1), vy(ind1),  dz*smart_deriv<1,0,0>(vx, ind1) + dz*smart_deriv<0,1,0>(vy, ind1) + vz.move(ind1, 0,0,1));
            bc.v_points.add_point(ind1, new_vec1);

            const auto ind2 = bc.old_new[bc.global_grid.convert_indices_unif(vec3(i,j,P) )];
            //const auto new_vec2 = vec3( -dx*smart_deriv<0,1,0>(vy, ind1) + -dx*smart_deriv<0,0,1>(vz, ind1) + vz.move(ind1, -1,0,0), vy(ind1), vz(ind1));
            const auto new_vec2 = vec3( vx(ind2), vy(ind2),  -dz*smart_deriv<1,0,0>(vx, ind2) - dz*smart_deriv<0,1,0>(vy, ind2) + vz.move(ind2, 0,0,-1));
            bc.v_points.add_point(ind2, new_vec2 );
        }
    }

    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            //need to change v.x
            const auto ind1 = bc.old_new[bc.global_grid.convert_indices_unif(vec3(0,j,k) )];
            //const auto new_vec1 = vec3( vx(ind1), dy*smart_deriv<1,0,0>(vx, ind1) + dy*smart_deriv<0,0,1>(vz, ind1) + vz.move(ind1, 0,1,0),  vz(ind1));
            const auto new_vec1 = vec3( dx*smart_deriv<0,1,0>(vy, ind1) + dx*smart_deriv<0,0,1>(vz, ind1) + vz.move(ind1, 1,0,0), vy(ind1), vz(ind1));
            bc.v_points.add_point( ind1, new_vec1 );

            const auto ind2 = bc.old_new[bc.global_grid.convert_indices_unif(vec3(N,j,k) )];
            //const auto new_vec2 = vec3( vx(ind1), -dy*smart_deriv<1,0,0>(vx, ind1) - dy*smart_deriv<0,0,1>(vz, ind1) + vz.move(ind1, 0,-1,0),  vz(ind1));
            const auto new_vec2 = vec3( -dx*smart_deriv<0,1,0>(vy, ind2) + -dx*smart_deriv<0,0,1>(vz, ind2) + vz.move(ind2, -1,0,0), vy(ind2), vz(ind2));
            bc.v_points.add_point(ind2, new_vec2);
        }
    }

    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            //need to change v.y
            const auto ind1 = bc.old_new[bc.global_grid.convert_indices_unif(vec3(i,0,k) )];
            //const auto new_vec1 = vec3( vx(ind1), vy(ind1),  dz*smart_deriv<1,0,0>(vx, ind1) + dz*smart_deriv<0,1,0>(vy, ind1) + vz.move(ind1, 0,0,1));
            const auto new_vec1 = vec3( vx(ind1), dy*smart_deriv<1,0,0>(vx, ind1) + dy*smart_deriv<0,0,1>(vz, ind1) + vz.move(ind1, 0,1,0),  vz(ind1));
            bc.v_points.add_point(ind1, new_vec1 );

            const auto ind2 = bc.old_new[bc.global_grid.convert_indices_unif(vec3(i,M,k) )];
            //const auto new_vec2 = vec3( vx(ind1), vy(ind1),  -dz*smart_deriv<1,0,0>(vx, ind1) - dz*smart_deriv<0,1,0>(vy, ind1) + vz.move(ind1, 01,0,-1));
            const auto new_vec2 = vec3( vx(ind2), -dy*smart_deriv<1,0,0>(vx, ind2) - dy*smart_deriv<0,0,1>(vz, ind2) + vz.move(ind2, 0,-1,0),  vz(ind2));
            bc.v_points.add_point(ind2, new_vec2 );
        }
    }

}*/

template <bool err>
bool update_pressure_BC(const boundary_conditions &BC, big_vec_d &p, const double accuracy_percent = 1.0) {
        for (unsigned i = 0; i < p.size(); i++) {
            if (p.is_boundary(i)) {
                const auto old = p(i);

                if (!BC.global_grid.off_walls(i)) {
                    p(i) = 0;
                } else {

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
                        } else {   //forward difference
                            p(i) = ny/nx* smart_deriv<0,1,0>(p, i)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i)*2*dx/3 - p.move(i,2,0,0)/3 + 4*p.move(i,1,0,0)/3;

                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p.has_up(i)) {  //backward difference
                            p(i) =-nx/ny* smart_deriv<1,0,0>(p, i)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i)*2*dy/3 - p.move(i,0,-2,0)/3 + 4*p.move(i,0,-1,0)/3;
                        } else { //forward difference
                            p(i) = nx/ny* smart_deriv<1,0,0>(p, i)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i)*2*dy/3 - p.move(i, 0,2,0)/3 + 4*p.move(i,0,1,0)/3;
                        }
                    } else {    //z direction
#ifndef NDEBUG
                        if (big_dir > 2) {
                            std::cerr << "the biggest direction cannot be larger than 2\n";
                        }
#endif
                        if (!p.has_back(i)) { //backward difference
                            p(i) =-ny/nz* smart_deriv<0,1,0>(p, i)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i)*2*dz/3 - p.move(i,0,0,-2)/3 + 4*p.move(i,0,0,-1)/3;
                        } else {  //forward difference
                            p(i) = ny/nz* smart_deriv<0,1,0>(p, i)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i)*2*dz/3 - p.move(i,0,0,2)/3 + 4*p.move(i,0,0,1)/3;

                        }

                    }


#ifdef UPDATE_VECS_CHECK_RESULTS_LOG
                    if (!std::isfinite(p(i))) {
                        std::cerr << "pressure boundary condition returned an infinite value\n";
                    }
                    if constexpr (err) {
                        if (std::abs(old-p(i))/std::abs(p(i))*100 > accuracy_percent) {
                            std::cerr << "more than 1% correction on pressure. This is bad\n";
                            if (big_dir == 0) {
                                std::cerr << "\tThis happened for x the largest direction\n";
                            } else if (big_dir == 1) {
                                std::cerr << "\tThis happened for y the largest direction\n";
                            } else {
                                std::cerr << "\tThis happened for z the largest direction\n";
                            }

                            std::cerr << "\tcorrection : " << std::abs(old-p(i))/std::abs(p(i))*100 << "%\n";
                            std::cerr << "\tat index " << i << "\n";
                            std::cerr << "\told pressure = " << old << "\t new pressure = " << p(i) << "\n";
                            std::cerr << "\tnormal vector = (" << norm << ")\n";

                            std::cerr << "\tThis was done using : ";
                            if (!p.has_right(i)) {
                                std::cerr << "backwards in x,  ";
                            } else {
                                std::cerr << "forwards in x,  ";
                            }
                            if (!p.has_up(i)) {
                                std::cerr << "backwards in y,  ";
                            } else {
                                std::cerr << "forwards in y,  ";
                            }
                            if (!p.has_back(i)) {
                                std::cerr << "backwards in z";
                            } else {
                                std::cerr << "forwards in z";
                            }
                            std::cerr << "\n";
                        }
                    }
#endif

                }
                if constexpr (err) {
                    if (std::abs(old-p(i))/std::abs(p(i))*100 > accuracy_percent) {
                        return false;
                    }
                }
            }

        }


    return true;
}

#endif //CODE_UPDATE_VECS_HPP
