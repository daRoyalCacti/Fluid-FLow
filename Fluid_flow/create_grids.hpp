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
    for (int k = 0; k < sz; k++) {
        for (int j = 0; j < sy; j++) {
            for (int i = 0; i < sx; i++) {


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


            const auto ind1 = g.convert_indices_unif(inds1);
            const auto ind2 = g.convert_indices_unif(inds2);

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
                inside_indices.insert( g.convert_indices_unif(ind_cpy) );
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

    g.create_no_points_unif();

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

    //removing points
    //======================================================
    std::vector<double> x, y, z;
    x.resize(g.x.size());
    y.resize(g.y.size());
    z.resize(g.z.size());
    std::vector<grid_relation> r;
    r.resize(g.r.size());
    std::map<unsigned, int> old_new;

    unsigned index_counter = 0;
    for (unsigned i = 0; i < g.x.size(); i++) {
        if (inside_indices.contains(i)) {
            old_new.insert({i, -1});
        } else {
            //std::cerr << index_counter << "\n";
            x[index_counter] = g.x[i];
            y[index_counter] = g.y[i];
            z[index_counter] = g.z[i];
            r[index_counter] = g.r[i];

            if (inside_indices.contains(r[index_counter].left)) { r[index_counter].left = -1; }
            if (inside_indices.contains(r[index_counter].right)) { r[index_counter].right = -1; }
            if (inside_indices.contains(r[index_counter].down)) { r[index_counter].down = -1; }
            if (inside_indices.contains(r[index_counter].up)) { r[index_counter].up = -1; }
            if (inside_indices.contains(r[index_counter].front)) { r[index_counter].front = -1; }
            if (inside_indices.contains(r[index_counter].back)) { r[index_counter].back -1; }
            old_new.insert({i, index_counter});
            index_counter++;
        }

    }
    x.shrink_to_fit();
    y.shrink_to_fit();
    z.shrink_to_fit();
    r.shrink_to_fit();
    g.x = std::move(x);
    g.y = std::move(y);
    g.z = std::move(z);
    g.r = std::move(r);

    boundary_normals norms_c(norms.size());
    std::cerr << "removing norms\n";
    for (const auto & n : norms.m) {
        norms_c.add_point( old_new.at(n.first), n.second  );
        //norms_c.add_point( old_new[n.first], n.second  );
    }
    norms = std::move(norms_c);

    mesh_points m_points_c(m_points.size());
    std::cerr << "removing points\n";
    for (const auto & n : m_points.m) {
        m_points_c.add_point( old_new.at(n.first), n.second  );
        //m_points_c.add_point( old_new[n.first], n.second  );
    }
    m_points = std::move(m_points_c);

    vel_points v_points_c(v_points.size());
    std::cerr << "removing vels\n";
    for (const auto & n : v_points.m) {
        v_points_c.add_point( old_new.at(n.first), n.second  );
        //v_points_c.add_point( old_new[n.first], n.second  );
    }
    v_points = std::move(v_points_c);

    //boundary_normals norms_c(norms.size());
    //mesh_points m_points_c(m_points.size());
    //vel_points v_points_c(v_points.size());

    //actually need to remove the points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

}



#endif //CODE_CREATE_GRIDS_HPP
