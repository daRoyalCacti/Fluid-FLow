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


//#define REMOVE_INSIDE_DLOG

void make_entire_grid(grid &g, const double Wx, const double Wy, const double Wz, const double dx, const double dy, const double dz, const double minx = 0, const double miny = 0, const double minz = 0) {
    const auto sx = static_cast<unsigned long>(Wx/dx);
    const auto sy = static_cast<unsigned long>(Wy/dy);
    const auto sz = static_cast<unsigned long>(Wz/dz);

    const vec3 mins = {minx,miny,minz};
    const vec3 maxs = {minx+Wx, miny+Wy, minz+Wz};
    /*
    g.dx = dx;
    g.dy = dy;
    g.dz = dz;*/
    g = grid(dx, dy, dz, mins, maxs);
    //std::cerr << g.dx << " " << g.dy << " " << g.dz << "\n";
    //std::cerr << dx << " " << dy << " " << dz << "\n";

    const auto s = static_cast<unsigned long>(sx*sy*sz);
    g.x.resize(s);
    g.y.resize(s);
    g.z.resize(s);
    g.r.resize(s);

    //g.create_no_points_unif();

    unsigned counter = 0;
    for (int k = 0; k < sz; k++) {
        for (int j = 0; j < sy; j++) {
            for (int i = 0; i < sx; i++) {


                g.x[counter] = minx + i*dx;
                g.y[counter] = miny + j*dy;
                g.z[counter] = minz + k*dz;

                int left, right, up, down, front, back;
                if (i != 0) {
                    left = static_cast<int>( g.convert_indices_unif(vec3(i-1,j,k)) );
                } else {
                    left = -1;
                }
                if (i != sx-1) {
                    right = static_cast<int>( g.convert_indices_unif(vec3(i+1,j,k)) );
                } else {
                    right = -1;
                }

                if (j != 0) {
                    down = static_cast<int>( g.convert_indices_unif(vec3(i,j-1,k)) );
                } else {
                    down = -1;
                }
                if (j != sy-1) {
                    up = static_cast<int>( g.convert_indices_unif(vec3(i,j+1,k)) );
                } else {
                    up = -1;
                }

                if (k != 0) {
                    front = static_cast<int>( g.convert_indices_unif(vec3(i,j,k-1)) );
                } else {
                    front = -1;
                }
                if (k != sz-1) {
                    back = static_cast<int>( g.convert_indices_unif(vec3(i,j,k+1)) );
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
                             boundary_normals &norms, mesh_points &m_points, vel_points &v_points, std::unordered_set<unsigned> &boundary_indices,
                             tri_inds& t_inds) noexcept {
    const auto hits = tm.get_collision_points(r);
    if (!hits.empty()) { //there was a collision
        for (auto h = hits.begin(); h != hits.end(); ++(++h)) {
            auto th = h;
            const auto col1 = th->second.v1;
            const auto norm1 = th->second.v2;
            const auto vel1 = th->second.v3;
            const triangle* t1 = th->second.v4;

            ++th;
            const auto col2 = th->second.v1;
            const auto norm2 = th->second.v2;
            const auto vel2 = th->second.v3;
            const triangle* t2 = th->second.v4;


            /*const auto inds1 = g.get_ind_unif(col1);
            const auto inds2 = g.get_ind_unif(col2);*/
            const auto inds1 = vec3{floor( (col1.x()-g.edge1.x()) /g.dx), floor( (col1.y()-g.edge1.y()) /g.dy), floor( (col1.z()-g.edge1.z())/g.dz)};
            const auto inds2 = vec3{floor( (col2.x()-g.edge1.x()) /g.dx), floor( (col2.y()-g.edge1.y()) /g.dy), floor( (col2.z()-g.edge1.z())/g.dz)};
            //std::cerr << inds1 << " " << inds2 << "\n";
            //std::cerr << g.dx << " " << col1.x() << " " << g.edge1 << "\n";


            const auto ind1 = g.convert_indices_unif(inds1);
            const auto ind2 = g.convert_indices_unif(inds2);


            unsigned axis;
            if (r.dir.x() == 1) {
                axis = 0;
#ifndef NDEBUG
                if (inds1.x() == inds2.x()) {
                    std::cerr << "ray hit the same point twice\n";
                }
#endif
            } if (r.dir.y() == 1) {
                axis = 1;
#ifndef NDEBUG
                if (inds1.y() == inds2.y()) {
                    std::cerr << "ray hit the same point twice\n";
                }
#endif
            } if (r.dir.z() == 1) {
                axis = 2;
#ifndef NDEBUG
                if (inds1.z() == inds2.z()) {
                    std::cerr << "ray hit the same point twice\n";
                }
#endif
            }



            //setting boundary points
            auto ind_cpy = inds1;
            for (auto i = static_cast<unsigned>(inds1[axis]); i <= inds2[axis]; i++) {
                ind_cpy[axis] = i;
                const auto insert_ind = g.convert_indices_unif(ind_cpy);
                inside_indices.insert( insert_ind );
            }

            //copying the data into the required data structures
            //===============================================================

            //if 2 rays collide with the mesh at the same point, take the average of the normals
            // - if more than 2 rays hits, the averaging weights the earlier hits weaker than the later hits (not desirable)
            // - doesn't seem to be working correctly
            if (norms.contains(ind1)) {
                norms.add_point(ind1, (norms.normal(ind1) + norm1)/2 );
            } else {
                norms.add_point(ind1, norm1);
            }

            if (norms.contains(ind2)) {
                norms.add_point(ind2, (norms.normal(ind2) + norm2)/2 );
            } else {
                norms.add_point(ind2, norm2);
            }



            //don't want to average collision point
            // - want to be sure it is on the mesh
            m_points.add_point(ind1, col1);
            m_points.add_point(ind2, col2);

            if (v_points.contains(ind1)) {
                v_points.add_point(ind1, (v_points.get_vel(ind1) + vel1)/2 );
            } else {
                v_points.add_point(ind1, vel1);
            }

            if (v_points.contains(ind2)) {
                v_points.add_point(ind2, (v_points.get_vel(ind2) + vel2)/2 );
            } else {
                v_points.add_point(ind2, vel2);
            }

            boundary_indices.insert(ind1);
            boundary_indices.insert(ind2);

            t_inds.add_point(ind1, t1);
            t_inds.add_point(ind2, t2);

        }

    }
}

//moves points inside boundary on a boring grid
void remove_inside_boundary_unif(grid &g, const triangle_mesh &tm, const mesh& m, boundary_normals &norms, mesh_points &m_points, vel_points &v_points, std::map<unsigned, int> &old_new, tri_inds& t_inds) {
    std::unordered_set<unsigned> inside_indices;
    std::unordered_set<unsigned> boundary_indices;

    constexpr unsigned test = 1;

    const auto minp = m.bounds.min;
    const auto maxp = m.bounds.max;

    const auto dims = ceil( (maxp - minp) / vec3(g.dx, g.dy, g.dz) );
    const auto Nx = static_cast<unsigned>(dims.x()) * test;
    const auto Ny = static_cast<unsigned>(dims.y()) * test;
    const auto Nz = static_cast<unsigned>(dims.z()) * test;


    const double dx = ( maxp.x() - minp.x() ) / ceil( ( maxp.x() - minp.x() ) /g.dx );
    const double dy = ( maxp.y() - minp.y() ) / ceil( ( maxp.y() - minp.y() ) /g.dy );
    const double dz = ( maxp.z() - minp.z() ) / ceil( ( maxp.z() - minp.z() ) /g.dz );

#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "rays from the z direction\n";
#endif
    //for (double x = minp.x(); x <= maxp.x(); x += g.dx/test) {
        //for (double y = minp.y(); y <= maxp.y(); y += g.dy/test) {
    for (unsigned i = 0; i <= Nx; i++) {
        for (unsigned j = 0; j <= Ny; j++) {
            //std::cerr << i << " " << j << "\n";
            const double x = minp.x() + i*dx/test;
            const double y = minp.y() + j*dy/test;

            const ray r(vec3(x, y, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points,boundary_indices, t_inds);
        }
    }

#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "rays from the y direction\n";
#endif
    //for (double x = minp.x(); x <= maxp.x(); x += g.dx/test) {
     //   for (double z = minp.z(); z <= maxp.z(); z += g.dz/test) {
     for (unsigned i = 0; i <= Nx; i++) {
         for (unsigned k = 0; k <= Nz; k++) {
             const double x = minp.x() + i*dx/test;
             const double z = minp.z() + k*dz/test;
            const ray r(vec3(x, 0, z), vec3(0,1,0) );
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points,boundary_indices, t_inds);
        }
    }

#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "rays from the x direction\n";
#endif
     for (unsigned j = 0; j <= Ny; j++) {
         for (unsigned k = 0; k <= Nz; k++) {
             const double y = minp.y() + j*dy/test;
             const double z = minp.z() + k*dz/test;
            const ray r(vec3(0, y, z), vec3(1,0,0) );
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points,boundary_indices, t_inds);
        }
    }

    //removing points
    //======================================================
#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "removing points\n";
#endif
    std::vector<double> x, y, z;
    x.reserve(g.x.size());
    y.reserve(g.y.size());
    z.reserve(g.z.size());
    std::vector<grid_relation> r;
    r.reserve(g.r.size());
    //std::map<unsigned, int> old_new;    //conversion between the old indices and the new (points removed) indices

    unsigned index_counter = 0;

    for (unsigned i = 0; i < g.x.size(); i++) {
        if (inside_indices.contains(i) && !boundary_indices.contains(i)) {
            old_new.insert({i, -1});
        } else {
            x.push_back(g.x[i]);
            y.push_back(g.y[i]);
            z.push_back(g.z[i]);
            r.push_back(g.r[i]);

            if (inside_indices.contains(r[index_counter].left) && !boundary_indices.contains(r[index_counter].left)) { r[index_counter].left = -1; }
            if (inside_indices.contains(r[index_counter].right) && !boundary_indices.contains(r[index_counter].right)) { r[index_counter].right = -1; }
            if (inside_indices.contains(r[index_counter].down) && !boundary_indices.contains(r[index_counter].down)) { r[index_counter].down = -1; }
            if (inside_indices.contains(r[index_counter].up) && !boundary_indices.contains(r[index_counter].up)) { r[index_counter].up = -1; }
            if (inside_indices.contains(r[index_counter].front) && !boundary_indices.contains(r[index_counter].front)) { r[index_counter].front = -1; }
            if (inside_indices.contains(r[index_counter].back) && !boundary_indices.contains(r[index_counter].back)) { r[index_counter].back = -1; }
            old_new.insert({i, index_counter});
            index_counter++;    //only gets incremented if not inside boundary
        }

    }
    x.shrink_to_fit();
    y.shrink_to_fit();
    z.shrink_to_fit();
    r.shrink_to_fit();

#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "recalibrating r\n";
#endif
    for (auto& e : r) {
        if (e.left != -1) { e.left = old_new.at(e.left); }
        if (e.right != -1) { e.right = old_new.at(e.right); }
        if (e.down != -1) { e.down = old_new.at(e.down); }
        if (e.up != -1) { e.up = old_new.at(e.up); }
        if (e.front != -1) {e.front = old_new.at(e.front);}
        if (e.back != -1) { e.back = old_new.at(e.back); }
    }


    g.x = std::move(x);
    g.y = std::move(y);
    g.z = std::move(z);
    g.r = std::move(r);


    boundary_normals norms_c(norms.size());
#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "recalibrating norms\n";
#endif
    for (const auto & n : norms.m) {

#ifndef NDEBUG
        norms_c.add_point( old_new.at(n.first), n.second  );
        if (old_new.at(n.first) == -1)  {
            std::cerr << "normal lies inside mesh\n";
        }
#else
        norms_c.add_point( old_new[n.first], n.second  );
#endif
    }
    norms = std::move(norms_c);

    mesh_points m_points_c(m_points.size());
#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "recalibrating points\n";
#endif
    for (const auto & n : m_points.m) {
#ifndef NDEBUG
        m_points_c.add_point( old_new.at(n.first), n.second  );
        if (old_new.at(n.first) == -1)  {
            std::cerr << "point lies inside mesh\n";
        }
#else
        m_points_c.add_point( old_new[n.first], n.second  );
#endif

    }
    m_points = std::move(m_points_c);

    vel_points v_points_c(v_points.size());
#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "recalibrating vels\n";
#endif
    for (const auto & n : v_points.m) {

#ifndef NDEBUG
        v_points_c.add_point( old_new.at(n.first), n.second  );
        if (old_new.at(n.first) == -1)  {
            std::cerr << "velocity lies inside mesh\n";
        }
#else
        v_points_c.add_point( old_new[n.first], n.second  );
#endif
    }
    v_points = std::move(v_points_c);

    tri_inds t_inds_c(t_inds.size());
#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "recalibrating triangle indices\n";
#endif
    for (const auto & n : t_inds.m) {
#ifndef NDEBUG
        t_inds_c.add_point( old_new.at(n.first), n.second  );
        if (old_new.at(n.first) == -1)  {
            std::cerr << "index lies inside mesh\n";
        }
#else
        t_inds_c.add_point( old_new[n.first], n.second  );
#endif
    }
    t_inds = std::move(t_inds_c);

    //std::cerr << t_inds.size() << " " << t_inds_c.size() << " " << v_points.size() << "\n";

#ifndef NDEBUG
#ifdef REMOVE_INSIDE_DLOG
        std::cerr << "checking results\n";
#endif
        //checking that all normals below to a boundary point
        /*for (const auto & n : norms.m) {
            if ( !g.is_boundary(n.first) ) {
                std::cerr << "normal at index " << n.first << " is not on a boundary point\n";
                std::cerr << "\tThis is at (" << g.x[n.first] << " " << g.y[n.first] << " " << g.z[n.first] << ")\n";
            }
        }*/

        //checking that all boundary points have a normal
        for (unsigned i = 0; i < g.r.size(); i++) {
            if (g.is_boundary(i) && g.off_walls(i) && !norms.contains(i)) {
                std::cerr << "boundary point at index " << i << " does not have a normal\n";
                std::cerr << "\tThis is at (" << g.x[i] << " " << g.y[i] << " " << g.z[i] << ")\n";
            }
        }

        /*for (const auto & n : v_points.m) {
            if ( !g.is_boundary(n.first) ) {
                std::cerr << "velocity at index " << n.first << " is not on a boundary point\n";
                std::cerr << "\tThis is at (" << g.x[n.first] << " " << g.y[n.first] << " " << g.z[n.first] << ")\n";
            }
        }*/

        //checking that all boundary points have a normal
        for (unsigned i = 0; i < g.r.size(); i++) {
            if (g.is_boundary(i) && g.off_walls(i) && !v_points.contains(i)) {
                std::cerr << "boundary point at index " << i << " does not have a velocity\n";
                std::cerr << "\tThis is at (" << g.x[i] << " " << g.y[i] << " " << g.z[i] << ")\n";
            }
        }

        //checking norms, vels, points, and inds all have the same points
        for (const auto &v : v_points.m) {
            if (!m_points.contains(v.first)) {
                std::cerr << "m_points does not contain " << v.first << " where v_points does\n";
            }
            if (!norms.contains(v.first)) {
                std::cerr << "norms does not contain " << v.first << " where v_points does\n";
            }
            if (!t_inds.contains(v.first)) {
                std::cerr << "t_inds does not contain " << v.first << " where v_points does\n";
            }
        }

        for (const auto &v : m_points.m) {
            if (!v_points.contains(v.first)) {
                std::cerr << "v_points does not contain " << v.first << " where m_points does\n";
            }
            if (!norms.contains(v.first)) {
                std::cerr << "norms does not contain " << v.first << " where m_points does\n";
            }
            if (!t_inds.contains(v.first)) {
                std::cerr << "t_inds does not contain " << v.first << " where m_points does\n";
            }
        }

        for (const auto &v : norms.m) {
            if (!v_points.contains(v.first)) {
                std::cerr << "v_points does not contain " << v.first << " where norms does\n";
            }
            if (!m_points.contains(v.first)) {
                std::cerr << "m_points does not contain " << v.first << " where norms does\n";
            }
            if (!t_inds.contains(v.first)) {
                std::cerr << "t_inds does not contain " << v.first << " where norms does\n";
            }
        }

        for (const auto &v : t_inds.m) {
            if (!m_points.contains(v.first)) {
                std::cerr << "m_points does not contain " << v.first << " where t_inds does\n";
            }
            if (!norms.contains(v.first)) {
                std::cerr << "norms does not contain " << v.first << " where t_inds does\n";
            }
            if (!v_points.contains(v.first)) {
                std::cerr << "v_points does not contain " << v.first << " where t_inds does\n";
            }
        }

#ifdef REMOVE_INSIDE_DLOG
        std::cerr << "finished checking results\n";
#endif
#endif


}



#endif //CODE_CREATE_GRIDS_HPP
