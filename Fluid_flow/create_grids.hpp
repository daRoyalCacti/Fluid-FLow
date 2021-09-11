//
// Created by jacob on 11/9/21.
//

#ifndef CODE_CREATE_GRIDS_HPP
#define CODE_CREATE_GRIDS_HPP

#include "../MyMath/grid.hpp"
#include <vector>
#include <unordered_set>

void make_enitre_gid(grid &g, const double Wx, const double Wy, const double Wz, const double dx, const double dy, const double dz) {
    const auto sx = static_cast<unsigned long>(Wx/dx);
    const auto sy = static_cast<unsigned long>(Wy/dy);
    const auto sz = static_cast<unsigned long>(Wz/dz);

    g.mins = {0,0,0};
    g.maxs = {Wx, Wy, Wz};

    const auto s = static_cast<unsigned long>(sx*sy*sz);
    g.x.resize(s);
    g.y.resize(s);
    g.z.resize(s);
    g.r.resize(s);

    unsigned counter = 0;
    for (int i = 0; i < sx; i++) {
        for (int j = 0; j < sy; j++) {
            for (int k = 0; k < sz; k++) {
                g.x[counter] = i*dx;
                g.y[counter] = j*dy;
                g.z[counter] = k*dz;

                int left, right, up, down, front, back;
                if (i != 0) {
                    left = i - 1;
                } else {
                    left = -1;
                }
                if (i != sx-1) {
                    right = i + 1;
                } else {
                    right = -1;
                }

                if (j != 0) {
                    down = j - 1;
                } else {
                    down = -1;
                }
                if (j != sy-1) {
                    up = j + 1;
                } else {
                    up = -1;
                }

                if (k != 0) {
                    front = k - 1;
                } else {
                    back = -1;
                }
                if (k != sz-1) {
                    front = k + 1;
                } else {
                    back = -1;
                }

                g.r[counter] = {left, right, down, up, front, back};

                counter++;
            }
        }
    }

}


void remove_inside_boundary(grid &g, const mesh &model, boundary_normals &norms, mesh_points &m_points) {
    std::unordered_set<unsigned> inside_indices;

    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned j = 0; j <= M; ++j) {
            const auto pos = vel_bc.get_pos(i,j,0);
            const auto dx = vel_bc.dx(i,j,0);
            const auto dy = vel_bc.dy(i,j,0);
            const ray r(vec3(pos.x() + dx/2, pos.y() + dy/2, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_z(r, inside_indices, norms, m_points, i, j);
        }
    }

    for (unsigned i = 0; i <= N; ++i) {
        for (unsigned k = 0; k <= P; ++k) {
            const auto pos = vel_bc.get_pos(i,0,k);
            const auto dx = vel_bc.dx(i,0,k);
            const auto dz = vel_bc.dz(i,0,k);
            const ray r(vec3(pos.x() + dx/2, 0, pos.z() + dz/2), vec3(0,1,0) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_y(r, inside_indices, norms, m_points, i, k);
        }
    }

    for (unsigned j = 0; j <= M; ++j) {
        for (unsigned k = 0; k <= P; ++k) {
            const auto pos = vel_bc.get_pos(0,j,k);
            const auto dy = vel_bc.dy(0,j,k);
            const auto dz = vel_bc.dz(0,j,k);
            const ray r(vec3(0, pos.y() + dy/2, pos.z() + dz/2), vec3(1,0,0) );//shoot ray through the middle of a grid point
            set_BC_mesh_1dir_x(r, inside_indices, norms, m_points, j, k);
        }
    }
}



#endif //CODE_CREATE_GRIDS_HPP
