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


    boundary_conditions() = delete;
    boundary_conditions(const mesh *m, const double dx, const double dy, const double dz) : tm(m) {
        std::cerr << "Remember to change velocity BC back to 0\n";

        set_wall_points();
        norms = boundary_normals<N,M,P>(num);
        create_wall_normals();

        p_bc = big_vec<N,M,P,double>(dx, dy, dz, &bound);
        vel_bc = big_vec<N,M,P,vec3>(dx, dy, dz, &bound);

        //needs to be called after p_bc is set because it uses dx, dy, dz from it
        update_mesh_boundary();

    }

    void update_pressure_BC();
    void update_velocity_BC();
    void update_mesh_boundary();

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

    void DEBUG_write_boundary_points() {
        std::ofstream output("../DEBUG/boundary_points.txt");
        if (output.is_open()) {
            for (unsigned i = 0; i <= N; i++) {
                for (unsigned j = 0; j <= M; j++) {
                    output << i * p_bc.dx << " " << j * p_bc.dy << " " << bound.is_boundary(i,j,P/2) << "\n";
                }
            }
        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }

    void DEBUG_write_normal_vectors() {
        std::ofstream output("../DEBUG/normal_vectors.txt");
        if (output.is_open()) {
            for (unsigned i = 0; i <= N; i++) {
                for (unsigned j = 0; j <= M; j++) {
                    output << i * p_bc.dx << " " << j * p_bc.dy << " ";
                    if (norms.contains(i,j,P/2)) {
                        output << norms.normal(i,j,P/2);
                    } else {
                        output << vec3(0);
                    }
                    output << "\n";
                }

            }
        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }

private:
    void set_wall_points();
    void create_wall_normals();
    void set_BC_mesh_1dir_z(const ray &r, std::vector<bool> &is_boundary, double dx, double dy, double dz, unsigned i, unsigned j) noexcept;
    void set_BC_mesh_1dir_y(const ray &r, std::vector<bool> &is_boundary, double dx, double dy, double dz, unsigned i, unsigned k) noexcept;
    void set_BC_mesh_1dir_x(const ray &r, std::vector<bool> &is_boundary, double dx, double dy, double dz, unsigned j, unsigned k) noexcept;
};


template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::set_BC_mesh_1dir_x(const ray &r, std::vector<bool> &is_boundary, const double dx, const double dy, const double dz, const unsigned j, const unsigned k) noexcept {
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

            const auto i_x1 = static_cast<unsigned>(floor(col1.x()/dx));
            const auto i_x2 = static_cast<unsigned>(floor(col2.x()/dx));

            //setting boundary points
            for (unsigned i = i_x1; i<=i_x2; i++) {
                is_boundary[i + (N+1)*j + (N+1)*(M+1)*k] = true;
            }

            //setting normals
            norms.add_point(i_x1, j, k, norm1);
            norms.add_point(i_x2, j, k, norm2);

            vel_bc.add_elm(i_x1, j, k, vel1);
            vel_bc.add_elm(i_x2, j, k, vel2);

        }

    }
}


template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::set_BC_mesh_1dir_y(const ray &r, std::vector<bool> &is_boundary, const double dx, const double dy, const double dz, const unsigned i, const unsigned k) noexcept {
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

            const auto i_y1 = static_cast<unsigned>(floor(col1.y()/dy));
            const auto i_y2 = static_cast<unsigned>(floor(col2.y()/dy));

            //setting boundary points
            for (unsigned j = i_y1; j<=i_y2; j++) {
                is_boundary[i + (N+1)*j + (N+1)*(M+1)*k] = true;
            }

            //setting normals
            norms.add_point(i, i_y1, k, norm1);
            norms.add_point(i, i_y2, k, norm2);

            vel_bc.add_elm(i, i_y1, k, vel1);
            vel_bc.add_elm(i, i_y2, k, vel2);

        }

    }
}


template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::set_BC_mesh_1dir_z(const ray &r, std::vector<bool> &is_boundary, const double dx, const double dy, const double dz, const unsigned i, const unsigned j) noexcept {
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

            const auto i_z1 = static_cast<unsigned>(floor(col1.z()/dz));
            const auto i_z2 = static_cast<unsigned>(floor(col2.z()/dz));

            //setting boundary points
            for (unsigned k = i_z1; k<=i_z2; k++) {
                is_boundary[i + (N+1)*j + (N+1)*(M+1)*k] = true;
            }

            //setting normals
            norms.add_point(i, j, i_z1, norm1);
            norms.add_point(i, j, i_z2, norm2);

            vel_bc.add_elm(i,j, i_z1, vel1);
            vel_bc.add_elm(i,j, i_z2, vel2);

        }

    }
}

template <unsigned N, unsigned M, unsigned P>
void boundary_conditions<N,M,P>::update_mesh_boundary() {
    //bool *is_boundary = new bool[(N+1)*(M+1)*(P+1)];
    //vector over c array so all elements initialised to 0 = false
    std::vector<bool> is_boundary;
    is_boundary.resize((N+1)*(M+1)*(P+1));

    //ray r(vec3{}, vec3(0,0,1));
    const auto dx = p_bc.dx;
    const auto dy = p_bc.dy;
    const auto dz = p_bc.dz;
    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned j = 0; j <= M; ++j) {
            const ray r(vec3(i*dx + dx/2, j*dy + dy/2, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_z(r, is_boundary, dx, dy, dz, i, j);

        }
    }

    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned k = 0; k <= P; ++k) {
            const ray r(vec3(i*dx + dx/2, 0, k*dz + dz/2), vec3(0,1,0) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_y(r, is_boundary, dx, dy, dz, i, k);
        }
    }

    for (unsigned j = 0; j <= M; ++j) {
        for (unsigned k = 0; k <= P; ++k) {
            const ray r(vec3(0, j*dy + dy/2, k*dz + dz/2), vec3(1,0,0) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_x(r, is_boundary, dx, dy, dz, j, k);
        }
    }

    //setting interior normals
    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned j = 0; j <= M; ++j) {
            for (unsigned k = 0; k <= P; ++k) {
                if (bound.is_boundary(i,j,k) && !norms.contains(i,j,k)) {
                    norms.add_point(i,j,k, vec3(0));
                }
            }
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


    //delete [] is_boundary;
}

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


#endif //CODE_BOUNDARY_CONDITIONS_HPP
