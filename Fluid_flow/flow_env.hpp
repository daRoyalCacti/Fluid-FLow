//
// Created by jacob on 21/8/21.
//

#ifndef CODE_FLOW_ENV_HPP
#define CODE_FLOW_ENV_HPP


template <unsigned N, unsigned M, unsigned P>
void set_BC(big_vec<N,M,P, vec3> &v, const double t) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            v.add_elm(i,j,0, 0,0,0);
            v.add_elm(i,j,P, 0,0,0);
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            v.add_elm(0,j,k, 0,0,0);
            v.add_elm(N,j,k, 0,0,0);
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            if (t < 0.2) {
                v.add_elm(i, 0, k, 0.1 * sin(0.1 * t), 0, 0);    //making the floor move
            } else {
                v.add_elm(i, 0, k, 0.1 * sin(0.1 * 0.2), 0, 0);
            }
            v.add_elm(i,M,k, 0,0,0);
        }
    }

}

template <unsigned N, unsigned M, unsigned P>
void enforce_PBC(big_vec<N,M,P, double> &p, const boundary_normals<N,M,P> &norms) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                if (p.is_boundary(i,j,k)) {

                    const auto norm = norms.normal(i,j,k);
                    const auto nx = norm.x();
                    const auto ny = norm.y();
                    const auto nz = norm.z();

                    const auto dx = p.dx;
                    const auto dy = p.dy;
                    const auto dz = p.dz;

                    //picking the direction
                    unsigned big_dir = 0;
                    if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                        big_dir = 1;
                    }
                    if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                        big_dir = 2;
                    }


                    if (big_dir == 0) { //x direction biggest
                        if (!p.has_left(i,j,k)) {   //forward difference
                            p(i,j,k) = ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 + nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i+2,j,k)/3 + 4*p(i+1,j,k)/3;
                        } else if (!p.has_right(i,j,k)) {   //backward difference
                            p(i,j,k) =-ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i-2,j,k)/3 + 4*p(i-1,j,k)/3;
                        } else {    //central difference
                            p(i+1,j,k) = -ny/nx * smart_deriv<0,1,0>(p, i,j,k)*dx - nz/nx * smart_deriv<0,0,1>(p, i,j,k)*dx + p(i-1,j,k);
                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p.has_down(i, j, k)) { //forward difference
                            p(i,j,k) = nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j+2,k)/3 + 4*p(i,j+1,k)/3;
                        } else if (!p.has_up(i, j, k)) {  //backward difference
                            p(i,j,k) =-nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j-2,k)/3 + 4*p(i,j-1,k)/3;
                        } else {  //central difference
                            p(i,j+1,k) = -nx/ny * smart_deriv<1,0,0>(p, i,j,k)*dy - nz/ny * smart_deriv<0,0,1>(p, i,j,k)*dy + p(i,j-1,k);
                        }
                    } else {    //z direction
#ifndef NDEBUG
                        if (big_dir > 2) {
                            std::cerr << "the biggest direction cannot be larger than 2\n";
                        }
#endif
                        if (!p.has_front(i,j,k)) {  //forward difference
                            p(i,j,k) = ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k+2)/3 + 4*p(i,j,k+1)/3;
                        } else if (!p.has_back(i,j,k)) { //backward difference
                            p(i,j,k) =-ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k-2)/3 + 4*p(i,j,k-1)/3;
                        } else {
                            p(i,j,k+1) = -ny/nz * smart_deriv<0,1,0>(p, i,j,k)*dz - nx/nz * smart_deriv<1,0,0>(p, i,j,k)*dz + p(i,j,k-1);
                        }
                    }

                }

            }
        }
    }
}

template <unsigned N, unsigned M, unsigned P>
void v_IC(big_vec<N,M,P, vec3> &v) {
    //just leave it as 0
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                v.add_elm(i, j, k, 0, 0, 0);
            }
        }

    }
}

//num is the number of boundary points
template <unsigned N, unsigned M, unsigned P>
void create_boundary_points(boundary_points<N,M,P> &bound, unsigned &num) {
    num = 0;
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
void create_boundary_normals(boundary_normals<N,M,P> &bound) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            bound.add_point(i,0,k, vec3(0,1,0));
            bound.add_point(i,M,k, vec3(0,-1,0));
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            bound.add_point(i,j,0, vec3(0,0,1));
            bound.add_point(i,j,P, vec3(0,0,-1));
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            bound.add_point(0,j,k,vec3(1,0,0));
            bound.add_point(N,j,k,vec3(-1,0,0));
        }
    }
}


template <unsigned N, unsigned M, unsigned P>
void DEBUG_check_normal_for_all_boundary_points(const boundary_points<N,M,P> &point, const boundary_normals<N,M,P> &norm) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <=M; j++) {
            for (unsigned k = 0; k<=P; k++) {
                if (point.is_boundary(i,j,k)) {
                    if (!norm.contains(i,j,k)) {
                        std::cerr << "Boundary at i=" << i << " j=" << j << " k=" << k << "does not have a normal\n";
                    }
                }
            }
        }
    }
}

#endif //CODE_FLOW_ENV_HPP
