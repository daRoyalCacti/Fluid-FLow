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

void get_mesh_collision(const triangle_mesh &tm, const grid &g, const ray &r, std::unordered_set<unsigned> &inside_indices, boundary_normals &norms, mesh_points &m_points) noexcept {
    const auto hits = tm.get_collision_points(r);
    if (!hits.empty()) { //there was a collision

        for (auto h = hits.begin(); h != hits.end(); ++(++h)) {
            auto th = h;
            const auto col1 = th->second.v1;
            const auto norm1 = th->second.v2;
            const auto vel1 = th->second.v3;

            ++th;
            const auto col2 = th->second.v1;
            const auto norm2 = th->second.v2;
            const auto vel2 = th->second.v3;

            const auto inds = g.get_inds(col1, col2);

            //const auto ind1 = g.get_ind(col1);
            //const auto ind2 = g.get_ind(col2);
            const auto ind1 = *inds.begin();
            const auto ind2 = *(--inds.end());

            //setting boundary points
            const auto begin_it = ++inds.begin();
            const auto end_it = --inds.end();
            for (auto i = begin_it; i != end_it; i++) {
                inside_indices.insert(*i);
            }

            //setting normals
            norms.add_point(ind1, norm1);
            norms.add_point(ind2, norm2);

            m_points.add_point(ind1, col1);
            m_points.add_point(ind2, col2);


        }

    }
}

void remove_inside_boundary(grid &g, const triangle_mesh &tm, boundary_normals &norms, mesh_points &m_points) {
    std::unordered_set<unsigned> inside_indices;

    for (double x = g.mins.x(); x < g.maxs.x(); x += g.dx) {
        for (double y = g.mins.y(); y < g.maxs.y(); y += g.dy) {
            const ray r(vec3(x+g.dx/2, y+g.dy/2, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            get_mesh_collision(tm, g, r, inside_indices, norms, m_points);
        }
    }
    for (double x = g.mins.x(); x < g.maxs.x(); x += g.dx) {
        for (double z = g.mins.z(); z < g.maxs.z(); z += g.dz) {
            const ray r(vec3(x+g.dx/2, 0, z + g.dz/2), vec3(0,1,0) );
            get_mesh_collision(tm, g, r, inside_indices, norms, m_points);
        }
    }
    for (double y = g.mins.y(); y < g.maxs.y(); y += g.dy) {
        for (double z = g.mins.z(); z < g.maxs.z(); z += g.dz) {
            const ray r(vec3(0, y+g.dy/2, z + g.dz/2), vec3(1,0,0) );
            get_mesh_collision(tm, g, r, inside_indices, norms, m_points);
        }
    }

}



#endif //CODE_CREATE_GRIDS_HPP
