//
// Created by jacob on 24/8/21.
//

#ifndef CODE_BOUNDARY_CONDITIONS_HPP
#define CODE_BOUNDARY_CONDITIONS_HPP


#include "../MyMath/boundary.hpp"
#include "../Rigid_body/body.hpp"
#include "../Rigid_body/triangle_mesh.hpp"
#include "../MyMath/grid.hpp"
#include "create_grids.hpp"
#include <cmath>


struct boundary_conditions {
    size_t no_inside_mesh = 0;
    grid global_grid;
    boundary_normals norms;
    mesh_points m_points;
    vel_points v_points;



    triangle_mesh tm;


    boundary_conditions() = delete;
    boundary_conditions(const mesh *m, const double dx, const double dy, const double dz, const double Wx, const double Wy, const double Wz, const double minx = 0, const double miny = 0, const double minz = 0) : tm(m) {
        std::cerr << "making entire grid\n";
        make_entire_grid(global_grid, Wx, Wy, Wz, dx, dy, dz, minx, miny, minz);
        //global_grid.create_no_points_unif();  //now called in make_entire grid

        std::cerr << "setting boundary normals\n";
        //set_wall_points();
        norms = boundary_normals(no_wall_points());
        std::cerr << "creating wall normals\n";
        create_wall_normals();
        std::cerr << "updating velocities at wall\n";
        update_velocity_wall_BC();

        //p_bc = big_vec<N,M,P,double>(dx, dy, dz, &bound);
        //vel_bc = big_vec<N,M,P,vec3>(dx, dy, dz, &bound);



        std::cerr << "removing inside points\n";
        remove_inside_boundary_unif(global_grid, tm, norms, m_points, v_points);

        //needs to be called after vel_bc is set because it uses dx, dy, dz from it
        /*
        update_mesh_boundary();
         */
    }

    //should never be used in real flow, only used for testing derivatives
    boundary_conditions(const double dx, const double dy, const double dz, const double Wx, const double Wy, const double Wz, const double minx = 0, const double miny = 0, const double minz = 0) : tm{} {
        std::cerr << "making entire grid\n";
        make_entire_grid(global_grid, Wx, Wy, Wz, dx, dy, dz, minx, miny, minz);
        //global_grid.create_no_points_unif();  //now called in make_entire grid

        std::cerr << "setting boundary normals\n";
        norms = boundary_normals(no_wall_points());
        std::cerr << "creating wall normals\n";
        create_wall_normals();
        std::cerr << "updating velocities at wall\n";
        update_velocity_wall_BC();
    }


    unsigned no_wall_points() const {
        const auto dims = global_grid.no_points_unif;
        const auto N = static_cast<unsigned>(dims.x()-1);
        const auto M = static_cast<unsigned>(dims.y()-1);
        const auto P = static_cast<unsigned>(dims.z()-1);
        return 2*(N+1)*(P+1) + 2*(N-1)*(M+1) + 2*(M-1)*(P-1);
    }
    /*

    void update_mesh_boundary();

    void extrapolate(big_vec<N,M,P, double> &p);
    void extrapolate(big_vec<N,M,P, vec3> &v) {
        extrapolate(v.xv);
        extrapolate(v.yv);
        extrapolate(v.zv);
    }

    //v is vector to enforce the conditions on
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
     */


    void DEBUG_write_boundary_points() const {
        std::ofstream output("../DEBUG/boundary_points.txt");
        if (output.is_open()) {
            const auto inds = global_grid.get_middle_inds();
            for (const auto ind : inds) {
                output << global_grid[ind].x() << " " << global_grid[ind].y() << " " << norms.contains(ind) << "\n";
            }

        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }


    void DEBUG_write_normal_vectors() const {
        std::ofstream output("../DEBUG/normal_vectors.txt");
        if (output.is_open()) {
            const auto inds = global_grid.get_middle_inds();
            for (const auto ind : inds) {
                output << global_grid[ind].x() << " " << global_grid[ind].y() << " ";
                if (norms.contains(ind)) {
                    output << norms.normal(ind);
                } else {
                    output << vec3(0);
                }
                output << "\n";
            }


        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }

    /*
    void update_pressure_BC(big_vec<N,M,P, double> &p);
     */
    void update_velocity_wall_BC();

private:
    void set_wall_points();
    void create_wall_normals();
    /*
    void set_BC_mesh_1dir_z(const ray &r, std::vector<bool> &is_boundary, unsigned i, unsigned j) noexcept;
    void set_BC_mesh_1dir_y(const ray &r, std::vector<bool> &is_boundary, unsigned i, unsigned k) noexcept;
    void set_BC_mesh_1dir_x(const ray &r, std::vector<bool> &is_boundary, unsigned j, unsigned k) noexcept;
    void set_matrix_row(unsigned x, unsigned y, unsigned z, unsigned &counter, Eigen::Matrix<double, 8, 8> &mat,  Eigen::Matrix<double, 8, 1> &vec, const big_vec<N,M,P, double> &p, double xi, double yi, double zi);
    */
};

/*
void boundary_conditions<N,M,P>::set_BC_mesh_1dir_x(const ray &r, std::vector<bool> &is_boundary, const unsigned  j, const unsigned k) noexcept {
    const auto hits = tm.get_collision_points(r);
    if (!hits.empty()) { //there was a collision

        for (auto h = hits.begin(); h != --hits.end(); ++h) {
            auto th = h;
            const auto col1 = th->second.v1;
            const auto norm1 = th->second.v2;
            const auto vel1 = th->second.v3;

            ++th;
            const auto col2 = th->second.v1;
            const auto norm2 = th->second.v2;
            const auto vel2 = th->second.v3;

            const auto ind1 = vel_bc.get_inds(col1);
            const auto ind2 = vel_bc.get_inds(col2);
            const auto i_x1 = ind1.x();
            const auto i_x2 = ind2.x();

            //setting boundary points
            for (unsigned i = i_x1; i<=i_x2; i++) {
                is_boundary[i + (N+1)*j + (N+1)*(M+1)*k] = true;
            }

            //setting normals
            norms.add_point(i_x1, j, k, norm1);
            norms.add_point(i_x2, j, k, norm2);

            vel_bc.add_elm(i_x1, j, k, vel1);
            vel_bc.add_elm(i_x2, j, k, vel2);

            points.add_point(i_x1, j, k, col1);
            points.add_point(i_x2, j, k, col2);


        }

    }
}


void boundary_conditions<N,M,P>::set_BC_mesh_1dir_y(const ray &r, std::vector<bool> &is_boundary, const unsigned i, const unsigned k) noexcept {
    const auto hits = tm.get_collision_points(r);
    if (!hits.empty()) { //there was a collision

        for (auto h = hits.begin(); h != --hits.end(); ++h) {
            auto th = h;
            const auto col1 = th->second.v1;
            const auto norm1 = th->second.v2;
            const auto vel1 = th->second.v3;

            ++th;
            const auto col2 = th->second.v1;
            const auto norm2 = th->second.v2;
            const auto vel2 = th->second.v3;

            const auto ind1 = vel_bc.get_inds(col1);
            const auto ind2 = vel_bc.get_inds(col2);
            const auto i_y1 = ind1.y();
            const auto i_y2 = ind2.y();

            //setting boundary points
            for (unsigned j = i_y1; j<=i_y2; j++) {
                is_boundary[i + (N+1)*j + (N+1)*(M+1)*k] = true;
            }

            //setting normals
            norms.add_point(i, i_y1, k, norm1);
            norms.add_point(i, i_y2, k, norm2);

            vel_bc.add_elm(i, i_y1, k, vel1);
            vel_bc.add_elm(i, i_y2, k, vel2);

            points.add_point(i, i_y1, k, col1);
            points.add_point(i, i_y2, k, col2);
        }

    }
}


void boundary_conditions<N,M,P>::set_BC_mesh_1dir_z(const ray &r, std::vector<bool> &is_boundary, const unsigned i, const unsigned j) noexcept {
    const auto hits = tm.get_collision_points(r);
    if (!hits.empty()) { //there was a collision

        for (auto h = hits.begin(); h != --hits.end(); ++h) {
            auto th = h;
            const auto col1 = th->second.v1;
            const auto norm1 = th->second.v2;
            const auto vel1 = th->second.v3;

            ++th;
            const auto col2 = th->second.v1;
            const auto norm2 = th->second.v2;
            const auto vel2 = th->second.v3;

            const auto ind1 = vel_bc.get_inds(col1);
            const auto ind2 = vel_bc.get_inds(col2);
            const auto i_z1 = ind1.z();
            const auto i_z2 = ind2.z();

            //setting boundary points
            for (unsigned k = i_z1; k<=i_z2; k++) {
                is_boundary[i + (N+1)*j + (N+1)*(M+1)*k] = true;
            }

            //setting normals
            norms.add_point(i, j, i_z1, norm1);
            norms.add_point(i, j, i_z2, norm2);

            vel_bc.add_elm(i,j, i_z1, vel1);
            vel_bc.add_elm(i,j, i_z2, vel2);

            points.add_point(i,j, i_z1, col1);
            points.add_point(i,j, i_z2, col2);
        }

    }
}


void boundary_conditions<N,M,P>::update_mesh_boundary() {
    //vector over c array so all elements initialised to 0 = false
    norms_prev = norms;
    bound.clear();
    set_wall_points();

    norms.clear();
    create_wall_normals();

    vel_bc.clear();
    update_velocity_wall_BC();
    points.clear();

    std::vector<bool> is_boundary;
    is_boundary.resize((N+1)*(M+1)*(P+1));

    //ray r(vec3{}, vec3(0,0,1));

    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned j = 0; j <= M; ++j) {
            const auto pos = vel_bc.get_pos(i,j,0);
            const auto dx = vel_bc.dx(i,j,0);
            const auto dy = vel_bc.dy(i,j,0);
            const ray r(vec3(pos.x() + dx/2, pos.y() + dy/2, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_z(r, is_boundary, i, j);
        }
    }

    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned k = 0; k <= P; ++k) {
            const auto pos = vel_bc.get_pos(i,0,k);
            const auto dx = vel_bc.dx(i,0,k);
            const auto dz = vel_bc.dz(i,0,k);
            const ray r(vec3(pos.x() + dx/2, 0, pos.z() + dz/2), vec3(0,1,0) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_y(r, is_boundary, i, k);
        }
    }

    for (unsigned j = 0; j <= M; ++j) {
        for (unsigned k = 0; k <= P; ++k) {
            const auto pos = vel_bc.get_pos(0,j,k);
            const auto dy = vel_bc.dy(0,j,k);
            const auto dz = vel_bc.dz(0,j,k);
            const ray r(vec3(0, pos.y() + dy/2, pos.z() + dz/2), vec3(1,0,0) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_x(r, is_boundary, j, k);
        }
    }

    //setting left and right and such
    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned j = 0; j <= M; ++j) {
            for (unsigned k = 0; k <= P; k++) {
                //note
                if (is_boundary[i + (N+1)*j + (N+1)*(M+1)*k]) { //dealing with boundary points
#ifndef NDEBUG
                    if (i <= 0) {
                        std::cerr << "setting BC from mesh: trying to access i < 0\n";
                    }
                    if (i >= N) {
                        std::cerr << "setting BC from mesh: trying to access i > N\n";
                    }
                    if (j <= 0) {
                        std::cerr << "setting BC from mesh: trying to access j < 0\n";
                    }
                    if (j >= M) {
                        std::cerr << "setting BC from mesh: trying to access j > M\n";
                    }
                    if (k <= 0) {
                        std::cerr << "setting BC from mesh: trying to access k < 0\n";
                    }
                    if (k >= P) {
                        std::cerr << "setting BC from mesh: trying to access k > P\n";
                    }
#endif
                    if (is_boundary[i+1 + (N+1)*j + (N+1)*(M+1)*k]) {
                        bound(i,j,k).has_left = false;
                    }
                    if (is_boundary[i-1 + (N+1)*j + (N+1)*(M+1)*k]) {
                        bound(i,j,k).has_right = false;
                    }
                    if (is_boundary[i + (N+1)*(j+1) + (N+1)*(M+1)*k]) {
                        bound(i,j,k).has_up = false;
                    }
                    if (is_boundary[i + (N+1)*(j-1) + (N+1)*(M+1)*k]) {
                        bound(i,j,k).has_down = false;
                    }
                    if (is_boundary[i + (N+1)*j + (N+1)*(M+1)*(k+1)]) {
                        bound(i,j,k).has_back = false;
                    }
                    if (is_boundary[i + (N+1)*j + (N+1)*(M+1)*(k-1)]) {
                        bound(i,j,k).has_front = false;
                    }




                }
            }
        }
    }

    no_inside_mesh = 0;
    //setting interior normals
    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned j = 0; j <= M; ++j) {
            for (unsigned k = 0; k <= P; ++k) {
                if (bound.is_boundary(i,j,k) && !norms.contains(i,j,k)) {
                    ++no_inside_mesh;
                    norms.add_point(i,j,k, vec3(0,0,0));
                }
            }
        }
    }

    //need to update the triangle mesh
    tm.update();


}

template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::set_wall_points() {
    //bottom and top
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(i,0,k).has_down = false;
            bound(i,M,k).has_up = false;
        }
    }
    //front and back
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            bound(i,j,0).has_front = false;
            bound(i,j,P).has_back = false;
        }
    }
    //top and bottom
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            bound(0,j,k).has_left = false;
            bound(N,j,k).has_right = false;
        }
    }
}



template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::update_pressure_BC(big_vec<N,M,P, double> &p) {
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                if (p.is_boundary(i,j,k)) {
#ifndef NDEBUG
                    if (!norms.contains(i,j,k)) {
                        std::cerr << "norms does not contain (" << i << ", " << j << ", " << k << ")\n";
                    }
#endif

                    const auto norm = norms.normal(i,j,k);

                    //this occurs when inside a boundary
                    if (norm == vec3(0) ) {
                        p(i,j,k) = 0;
                        continue;   //rest of the code will only error
                    }

                    const auto nx = norm.x();
                    const auto ny = norm.y();
                    const auto nz = norm.z();

                    const auto dx = p.dx(i,j,k);
                    const auto dy = p.dy(i,j,k);
                    const auto dz = p.dz(i,j,k);

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
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        } else if (!p.has_right(i,j,k)) {   //backward difference
                            p(i,j,k) =-ny/nx* smart_deriv<0,1,0>(p, i,j,k)*2*dx/3 - nz/nx* smart_deriv<0,0,1>(p, i,j,k)*2*dx/3 - p(i-2,j,k)/3 + 4*p(i-1,j,k)/3;
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        } else {    //central difference
                            p(i+1,j,k) = -ny/nx * smart_deriv<0,1,0>(p, i,j,k)*dx - nz/nx * smart_deriv<0,0,1>(p, i,j,k)*dx + p(i-1,j,k);
#ifndef NDEBUG
                            if (!std::isfinite(p(i+1,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        }
                    } else if (big_dir == 1) {  //y direction biggest
                        if (!p.has_down(i, j, k)) { //forward difference
                            p(i,j,k) = nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 + nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j+2,k)/3 + 4*p(i,j+1,k)/3;
 #ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        } else if (!p.has_up(i, j, k)) {  //backward difference
                            p(i,j,k) =-nx/ny* smart_deriv<1,0,0>(p, i,j,k)*2*dy/3 - nz/ny* smart_deriv<0,0,1>(p, i,j,k)*2*dy/3 - p(i,j-2,k)/3 + 4*p(i,j-1,k)/3;
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        } else {  //central difference
                            p(i,j+1,k) = -nx/ny * smart_deriv<1,0,0>(p, i,j,k)*dy - nz/ny * smart_deriv<0,0,1>(p, i,j,k)*dy + p(i,j-1,k);
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j+1,k))) {
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
                    if (!p.has_front(i,j,k)) {  //forward difference
                            p(i,j,k) = ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 + nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k+2)/3 + 4*p(i,j,k+1)/3;
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        } else if (!p.has_back(i,j,k)) { //backward difference
                            p(i,j,k) =-ny/nz* smart_deriv<0,1,0>(p, i,j,k)*2*dz/3 - nx/nz* smart_deriv<1,0,0>(p, i,j,k)*2*dz/3 - p(i,j,k-2)/3 + 4*p(i,j,k-1)/3;
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        } else {
                            p(i,j,k+1) = -ny/nz * smart_deriv<0,1,0>(p, i,j,k)*dz - nx/nz * smart_deriv<1,0,0>(p, i,j,k)*dz + p(i,j,k-1);
#ifndef NDEBUG
                            if (!std::isfinite(p(i,j,k+1))) {
                                std::cerr << "pressure boundary condition returned an infinite value\n";
                            }
#endif
                        }


                    }

                }

            }
        }
    }
}
*/

void boundary_conditions::update_velocity_wall_BC()  {
    const auto dims = global_grid.no_points_unif;
    const auto N = static_cast<unsigned>(dims.x()-1);
    const auto M = static_cast<unsigned>(dims.y()-1);
    const auto P = static_cast<unsigned>(dims.z()-1);

    //setting the walls
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            v_points.add_point(global_grid.convert_indices_unif(vec3(i,j,0) ), vec3(0) );
            v_points.add_point(global_grid.convert_indices_unif(vec3(i,j,P) ), vec3(0) );
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            v_points.add_point(global_grid.convert_indices_unif(vec3(0,j,k) ), vec3(0) );
            v_points.add_point(global_grid.convert_indices_unif(vec3(N,j,k) ), vec3(0) );
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
#ifdef MOVING_WALL
            v_points.add_point(global_grid.convert_indices_unif(vec3(i,0,k) ), vec3(1,0,0) );
#else
            v_points.add_point(global_grid.convert_indices_unif(vec3(i,0,k) ), vec3(0) );
            vel_bc.add_elm(i, 0, k, 0, 0, 0);
#endif
            v_points.add_point(global_grid.convert_indices_unif(vec3(i,M,k) ), vec3(0) );
        }
    }

}


void boundary_conditions::create_wall_normals() {
    const auto dims = global_grid.no_points_unif;
    const auto N = static_cast<unsigned>(dims.x()-1);
    const auto M = static_cast<unsigned>(dims.y()-1);
    const auto P = static_cast<unsigned>(dims.z()-1);

    for (unsigned i = 0; i <= N; i++) {
        for (unsigned k = 0; k <= P; k++) {
            norms.add_point(global_grid.convert_indices_unif(vec3(i,0,k) ), vec3(0,1,0));
            norms.add_point(global_grid.convert_indices_unif(vec3(i,M,k) ), vec3(0,-1,0));
        }
    }
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            norms.add_point(global_grid.convert_indices_unif(vec3(i,j,0) ), vec3(0,0,1));
            norms.add_point(global_grid.convert_indices_unif(vec3(i,j,P) ), vec3(0,0,-1));
        }
    }
    for (unsigned j = 0; j <= M; j++) {
        for (unsigned k = 0; k <= P; k++) {
            norms.add_point(global_grid.convert_indices_unif(vec3(0,j,k) ),vec3(1,0,0));
            norms.add_point(global_grid.convert_indices_unif(vec3(N,j,k) ),vec3(-1,0,0));
        }
    }
}
/*
//https://en.wikipedia.org/wiki/Trilinear_interpolation#Alternative_algorithm
// uses linear interpolation to update points that were previously inside a boundary
// would be quite simple to update to quadratic interpolation
// uses interpolation instead of an explicit extrapolation procedure is because I can't find any extrapolation procedures
template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::extrapolate(big_vec<N,M,P, double> &p) {
    //finds the constants in
    //y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz

    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <=P; k++) {
                //const bool is_bound = bound.is_boundary(i,j,k);
                //const bool was_bound = norms.contains(i,j,k) && norms.normal(i,j,k) != vec3(0); //bound_prev.is_boundary(i,j,k);

                const bool is_bound = norms.contains(i,j,k);
                const bool was_bound = norms_prev.contains(i,j,k) && norms_prev.normal(i,j,k) == vec3(0);


                if (!is_bound && was_bound ) {   //was inside mesh but is current outside of the boundary
                    //the points where the extrapolation is to be used
                    const unsigned xi = i;
                    const unsigned yi = j;
                    const unsigned zi = k;





                    Eigen::Matrix<double, 8, 8> mat;
                    Eigen::Matrix<double, 8, 1> vec;

                    //using values on the faces of ever-increasing cubes
                    for (int s = 2; s < 80; s+=2) {    //upper bound on s arbitrarily large
                        std::vector<unsigned> vals_b, vals_s;
                        vals_b.resize(s+1);
                        vals_s.resize(s-1);
                        std::iota(vals_b.begin(), vals_b.end(), -s/2);  //include the entire surface
                        std::iota(vals_s.begin(), vals_s.end(), -s/2+1);//surface less 1 edge

                        unsigned counter = 0;

                        //moving along the surfaces of the cubes
                        {
                            int x = -s/2;
                            for (const auto &y : vals_b) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int x = s/2;
                            for (const auto &y : vals_b) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int y = -s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int y = s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &z : vals_b) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int z = -s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &y : vals_s) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }
                        {
                            int z = s/2;
                            for (const auto &x : vals_s) {
                                for (const auto &y : vals_s) {
                                    set_matrix_row(x, y, z, counter, mat, vec, p, xi, yi, zi);
                                }
                            }
                        }

                        if (counter >= 8) {
                            break;
                        }


                    }   //end s for loop

                    //mat and vec now set, just need to solve for the coefficients
                    const Eigen::Matrix<double, 8, 1> a = mat.inverse()*vec;

                    //setting the value of p
                    //y = a0 + a1x+ a2y + a3z + a4xy + a5xz + a6yz + a7xyz
                    const auto pos = vel_bc.get_pos(xi,yi,zi);
                    const auto x = pos.x();
                    const auto y = pos.y();
                    const auto z = pos.z();
                    p(xi, yi, zi) = a(0) + a(1)*x + a(2)*y + a(3)*z + a(4)*x*y + a(5)*x*z + a(6)*y*z + a(7)*x*y*z;

                }

            }
        }
    }

}

template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::set_matrix_row(const unsigned x, const unsigned y, const unsigned z, unsigned &counter, Eigen::Matrix<double, 8, 8> &mat,
                                                Eigen::Matrix<double, 8, 1> &vec, const big_vec<N,M,P, double> &p, const double xi, const double yi, const double zi) {
    const auto xp = x+xi;
    const auto yp = y+yi;
    const auto zp = z+zi;

    std::cerr << "set matrix row will soon be wrong\n";


    //https://en.wikipedia.org/wiki/Trilinear_interpolation#Alternative_algorithm
    bool has_normal = norms.contains(xp, yp, zp);
    if (!has_normal || (norms.normal(xp, yp, zp) != vec3(0))  ) {    //if point not inside a boundary
        const auto pos = vel_bc.get_pos(xp, yp, zp);
        const auto x0 = pos.x();
        const auto y0 = pos.y();
        const auto z0 = pos.z();

        mat(counter, 0) = 1;
        mat(counter, 1) = x0;
        mat(counter, 2) = y0;
        mat(counter, 3) = z0;
        mat(counter, 4) = x0*y0;
        mat(counter, 5) = x0*z0;
        mat(counter, 6) = y0*z0;
        mat(counter, 7) = x0*y0*z0;

        vec(counter) = p(xp, yp, zp);
    }

    counter++;
}

*/

#endif //CODE_BOUNDARY_CONDITIONS_HPP
