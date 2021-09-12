//
// Created by jacob on 11/9/21.
//

#ifndef CODE_CREATE_GRIDS_HPP
#define CODE_CREATE_GRIDS_HPP

#include "../MyMath/grid.hpp"
#include "../Rigid_body/triangle_mesh.hpp"
#include "../Rigid_body/ray.hpp"
#include "../MyMath/boundary.hpp"
#include <vector>
#include <unordered_set>

void make_entire_grid(grid &g, const double Wx, const double Wy, const double Wz, const double dx, const double dy, const double dz) {
    const auto sx = static_cast<unsigned long>(Wx/dx);
    const auto sy = static_cast<unsigned long>(Wy/dy);
    const auto sz = static_cast<unsigned long>(Wz/dz);

    g.mins = {0,0,0};
    g.maxs = {Wx, Wy, Wz};
    g.dx = dx;
    g.dy = dy;
    g.dz = dz;

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

unsigned long convert_indices_unif(const vec3 &inds, const vec3 &no_points) {
    return static_cast<unsigned long>( inds.x() + inds.y()*no_points.x() + inds.z()*no_points.x()*no_points.y() );
}

void get_mesh_collision_unif(const triangle_mesh &tm, const grid &g, const ray &r, std::unordered_set<unsigned> &inside_indices,
                             boundary_normals &norms, mesh_points &m_points, vel_points &v_points) noexcept {
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


            //const auto ind1 = g.get_ind(col1);
            //const auto ind2 = g.get_ind(col2);
            //const auto ind1 = *inds.begin();
            //const auto ind2 = *(--inds.end());

            const auto inds1 = vec3( floor(col1.x()/g.dx), floor(col1.y()/g.dy),floor(col1.z()/g.dz) );
            const auto inds2 = vec3( floor(col2.x()/g.dx), floor(col2.y()/g.dy),floor(col2.z()/g.dz) );

            const vec3 no_points = round( (g.maxs - g.mins) / vec3(g.dx, g.dy, g.dz) );

            const auto ind1 = convert_indices_unif(inds1, no_points);
            const auto ind2 = convert_indices_unif(inds2, no_points);

            unsigned axis;
            if (r.dir.x() == 1) {
                axis = 0;
            } if (r.dir.y() == 1) {
                axis = 1;
            } if (r.dir.z() == 1) {
                axis = 2;
            }

            //setting boundary points
            auto ind_cpy = inds1;
             for (unsigned i = static_cast<unsigned>(inds1[axis]) + 1; i < inds2[axis]; i++) {
                ind_cpy[axis] = i;
                inside_indices.insert( convert_indices_unif(ind_cpy, no_points) );
            }

            //setting normals
            norms.add_point(ind1, norm1);
            norms.add_point(ind2, norm2);

            m_points.add_point(ind1, col1);
            m_points.add_point(ind2, col2);

            v_points.add_point(ind1, vel1);
            v_points.add_point(ind2, vel2);


        }

    }
}

//moves points inside boundary on a boring grid
void remove_inside_boundary_unif(grid &g, const triangle_mesh &tm, boundary_normals &norms, mesh_points &m_points, vel_points &v_points) {
    std::unordered_set<unsigned> inside_indices;

    for (double x = g.mins.x(); x < g.maxs.x(); x += g.dx) {
        for (double y = g.mins.y(); y < g.maxs.y(); y += g.dy) {
            const ray r(vec3(x+g.dx/2, y+g.dy/2, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points);
        }
    }
    for (double x = g.mins.x(); x < g.maxs.x(); x += g.dx) {
        for (double z = g.mins.z(); z < g.maxs.z(); z += g.dz) {
            const ray r(vec3(x+g.dx/2, 0, z + g.dz/2), vec3(0,1,0) );
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points);
        }
    }
    for (double y = g.mins.y(); y < g.maxs.y(); y += g.dy) {
        for (double z = g.mins.z(); z < g.maxs.z(); z += g.dz) {
            const ray r(vec3(0, y+g.dy/2, z + g.dz/2), vec3(1,0,0) );
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points);
        }
    }

}



#endif //CODE_CREATE_GRIDS_HPP
