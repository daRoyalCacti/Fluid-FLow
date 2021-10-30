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


#define REMOVE_INSIDE_DLOG
#define CHECK_GRID_RESULTS

#ifndef NDEBUG
#define CHECK_GRID_RESULTS
#endif

void make_entire_grid(grid &g, const double Wx, const double Wy, const double Wz, const double dx, const double dy, const double dz, const unsigned sx, const unsigned sy, unsigned sz, const double minx = 0, const double miny = 0, const double minz = 0) {
    /*const auto sx = static_cast<unsigned long>(Wx/dx);
    const auto sy = static_cast<unsigned long>(Wy/dy);
    const auto sz = static_cast<unsigned long>(Wz/dz);*/

    const vec3 mins = {minx,miny,minz};
    const vec3 maxs = {minx+Wx, miny+Wy, minz+Wz};

    g = grid(dx, dy, dz, mins, maxs, sx, sy, sz);

    const auto s = sx*sy*sz;
    g.x.resize(s);
    g.y.resize(s);
    g.z.resize(s);
    g.r.resize(s);


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
        if (hits.size() % 2 != 0) {
            //std::cerr << "Ray hit mesh an odd number of times\n";
            //std::cerr << "\t" << hits.size() << "\n";
            return;
        }
        for (auto h = hits.begin(); h != hits.end(); ++(++h)) {
            auto th = h;
            /*const auto col1 = th->second.v1;
            const auto norm1 = th->second.v2;
            const auto vel1 = th->second.v3;*/
            const triangle* t1 = th->second;
            const auto col1 = r.at(th->first);
            const auto norm1 = t1->get_normal(col1);
            const auto vel1 = t1->get_velocity(col1);

            /*if ( std::abs(norm1.length() - 1) < 0.001) {
                std::cerr << "norm1 is not unit\n";
            }*/

            ++th;
            /*const auto col2 = th->second.v1;
            const auto norm2 = th->second.v2;
            const auto vel2 = th->second.v3;*/
            const triangle* t2 = th->second;
            const auto col2 = r.at(th->first);
            const auto norm2 = t2->get_normal(col1);
            const auto vel2 = t2->get_velocity(col1);

            /*if ( std::abs(norm2.length() - 1) < 0.001) {
                std::cerr << "norm2 is not unit\n";
            }*/


            /*if (col1.z() == 1 || col2.z() == 1) {
                std::cerr << "\tthe hit\n";
            }*/


            //edge 1 is lower left corner
            const auto inds1 = vec3{floor( (col1.x()-g.edge1.x()) /g.dx), floor( (col1.y()-g.edge1.y()) /g.dy), floor( (col1.z()-g.edge1.z())/g.dz)};
            const auto inds2 = vec3{floor( (col2.x()-g.edge1.x()) /g.dx), floor( (col2.y()-g.edge1.y()) /g.dy), floor( (col2.z()-g.edge1.z())/g.dz)};


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
            //auto ind_cpy = inds1;
            inside_indices.insert( ind1 );
            auto temp_ind = ind1;
            while (temp_ind < ind2) {
                temp_ind = g.get_move_ind(temp_ind, r.dir);
                inside_indices.insert( temp_ind );
            }
            /*for (auto i = static_cast<unsigned>(inds1[axis]); i <= inds2[axis]; i++) {
                ind_cpy[axis] = i;
                const auto insert_ind = g.convert_indices_unif(ind_cpy);
                inside_indices.insert( insert_ind );
            }*/

            //copying the data into the required data structures
            //===============================================================

            //if 2 rays collide with the mesh at the same point, take the average of the normals
            // - if more than 2 rays hits, the averaging weights the earlier hits weaker than the later hits (not desirable)
            // - doesn't seem to be working correctly
            /*if (norms.contains(ind1)) {
                norms.add_point(ind1, (norms.normal(ind1) + norm1)/2 );
            } else {*/
                norms.add_point(ind1, norm1);
            //}

            /*if (norms.contains(ind2)) {
                norms.add_point(ind2, (norms.normal(ind2) + norm2)/2 );
            } else {*/
                norms.add_point(ind2, norm2);
            //}



            //don't want to average collision point
            // - want to be sure it is on the mesh
            m_points.add_point(ind1, col1);
            m_points.add_point(ind2, col2);

             /*if (v_points.contains(ind1)) {
                v_points.add_point(ind1, (v_points.get_vel(ind1) + vel1)/2 );
            } else {*/
                v_points.add_point(ind1, vel1);
            //}

            /*if (v_points.contains(ind2)) {
                v_points.add_point(ind2, (v_points.get_vel(ind2) + vel2)/2 );
            } else {*/
                v_points.add_point(ind2, vel2);
            //}

            boundary_indices.insert(ind1);
            boundary_indices.insert(ind2);

            t_inds.add_point(ind1, t1);
            t_inds.add_point(ind2, t2);

        }

    }
}


/*void update_r_point(grid_relation &r, const std::unordered_set<unsigned> &boundary_indices, const std::unordered_set<unsigned> &inside_indices) {
    if (inside_indices.contains(r.left) && !boundary_indices.contains(r.left)) { r.left = -1; }
    if (inside_indices.contains(r.right) && !boundary_indices.contains(r.right)) { r.right = -1; }
    if (inside_indices.contains(r.down) && !boundary_indices.contains(r.down)) { r.down = -1; }
    if (inside_indices.contains(r.up) && !boundary_indices.contains(r.up)) { r.up = -1; }
    if (inside_indices.contains(r.front) && !boundary_indices.contains(r.front)) { r.front = -1; }
    if (inside_indices.contains(r.back) && !boundary_indices.contains(r.back)) { r.back = -1; }

    if (boundary_indices.contains(r.left)) { r.left = -1; }
    if (boundary_indices.contains(r.right)) { r.right = -1; }
    if (boundary_indices.contains(r.down)) { r.down = -1; }
    if (boundary_indices.contains(r.up)) { r.up = -1; }
    if (boundary_indices.contains(r.front)) { r.front = -1; }
    if (boundary_indices.contains(r.back)) { r.back = -1; }
}*/

//moves points inside boundary on a boring grid
void remove_inside_boundary_unif(grid &g, const triangle_mesh &tm, const mesh& m, boundary_normals &norms, mesh_points &m_points, vel_points &v_points, std::map<unsigned, int> &old_new, tri_inds& t_inds) {
    std::unordered_set<unsigned> inside_indices;
    std::unordered_set<unsigned> boundary_indices;
    //auto & boundary_indices = g.boundary_indices;

    constexpr unsigned test = 4;

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
    for (unsigned i = 0; i <= Nx; i++) {
        for (unsigned j = 0; j <= Ny; j++) {
            const double x = minp.x() + i*dx/test + dx/2;
            const double y = minp.y() + j*dy/test + dy/2;

            const ray r(vec3(x, y, 0), vec3(0,0,1) );//shoot ray through the middle of a grid point
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points,boundary_indices, t_inds);
        }
    }

#ifdef REMOVE_INSIDE_DLOG
std::cerr << "rays from the y direction\n";
#endif
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
            const double y = minp.y() + j*dy/test + dy/2;
            const double z = minp.z() + k*dz/test + dz/2;
            const ray r(vec3(0, y, z), vec3(1,0,0) );
            get_mesh_collision_unif(tm, g, r, inside_indices, norms, m_points, v_points,boundary_indices, t_inds);
        }
    }



    //removing points
    //======================================================
#ifdef REMOVE_INSIDE_DLOG
std::cerr << "making derivatives work\n";
#endif
    std::vector<double> x, y, z;
    x.reserve(g.x.size());
    y.reserve(g.y.size());
    z.reserve(g.z.size());
    std::vector<grid_relation> r;
    r.reserve(g.r.size());




    //making singleton points boundary points (singleton as in you can't take a derivative there --- they are isolted)
    // not this is not guaranteed to fix the problem - making some points boundary could then mean that move points becomes singleton
    //  - try interating through 30 times to fix this
    bool was_update = true;
    while (was_update) {
        was_update = false;
        //unsigned counter = 0;
        //for (const auto i : boundary_indices) {
        //std::cerr << i << "\n";
        //std::cerr << counter++ << "/" << boundary_indices.size() << "\n";
        for (unsigned i = 0; i < g.size(); i++) {
            if (boundary_indices.contains(i)) {
                //std::cerr << i << "/" << g.size() << "\n";
                /*if (i==1023283) {
                    std::cerr < < g.r[i].left << "\n";
                    std::cerr << g.r[i].right << "\n";
                    std::cerr << g.r[i].down << "\n";
                    std::cerr << g.r[i].up << "\n";
                    std::cerr << g.r[i].front << "\n";
                    std::cerr << g.r[i].back << "\n";
                }*/


                //grid_relation r2 = g.r[i];
                //std::cerr << "\t1\n";
                //update_r_point(r2, boundary_indices, inside_indices);


                const bool will_left = (g.r[i].left != -1) && ( !inside_indices.contains(g.r[i].left)   || boundary_indices.contains(g.r[i].left));
                const bool will_right = (g.r[i].right != -1) && ( !inside_indices.contains(g.r[i].right)   || boundary_indices.contains(g.r[i].right));
                const bool will_down = (g.r[i].down != -1) && ( !inside_indices.contains(g.r[i].down)   || boundary_indices.contains(g.r[i].down));
                const bool will_up = (g.r[i].up != -1) && ( !inside_indices.contains(g.r[i].up)   || boundary_indices.contains(g.r[i].up));
                const bool will_front = (g.r[i].front != -1) && ( !inside_indices.contains(g.r[i].front)   || boundary_indices.contains(g.r[i].front));
                const bool will_back = (g.r[i].back != -1) && ( !inside_indices.contains(g.r[i].back)   || boundary_indices.contains(g.r[i].back));


                //std::cerr << "\t2\n";
                bool will_2left, will_2right, will_2down, will_2up, will_2front, will_2back;
                if (will_left) {
                    will_2left = (g.r[g.r[i].left].left != -1) && ( !inside_indices.contains(g.r[g.r[i].left].left) || boundary_indices.contains(g.r[g.r[i].left].left) );
                } else {
                    will_2left = false;
                }

                if (will_right) {
                    will_2right = (g.r[g.r[i].right].right != -1) && ( !inside_indices.contains(g.r[g.r[i].right].right) || boundary_indices.contains(g.r[g.r[i].right].right) );
                }else {
                    will_2right = false;
                }

                if (will_down) {
                    will_2down = (g.r[g.r[i].down].down != -1) && ( !inside_indices.contains(g.r[g.r[i].down].down) || boundary_indices.contains(g.r[g.r[i].down].down) );
                }else {
                    will_2down = false;
                }

                if (will_up) {
                    will_2up = (g.r[g.r[i].up].up != -1) && ( !inside_indices.contains(g.r[g.r[i].up].up) || boundary_indices.contains(g.r[g.r[i].up].up) );
                }else {
                    will_2up = false;
                }

                if (will_front) {
                    will_2front = (g.r[g.r[i].front].front != -1) && ( !inside_indices.contains(g.r[g.r[i].front].front) || boundary_indices.contains(g.r[g.r[i].front].front) );
                }else {
                    will_2front = false;
                }

                if (will_back) {
                    will_2back = (g.r[g.r[i].back].back != -1) && ( !inside_indices.contains(g.r[g.r[i].back].back) || boundary_indices.contains(g.r[g.r[i].back].back) );
                }else {
                    will_2back = false;
                }


                bool will_3left, will_3right, will_3down, will_3up, will_3front, will_3back;
                if (will_2left) {
                    will_3left = (g.r[g.r[g.r[i].left].left].left != -1) && ( !inside_indices.contains(g.r[g.r[g.r[i].left].left].left) || boundary_indices.contains(g.r[g.r[g.r[i].left].left].left) );
                } else {
                    will_3left = false;
                }

                if (will_2right) {
                    will_3right = (g.r[g.r[g.r[i].right].right].right != -1) && ( !inside_indices.contains(g.r[g.r[g.r[i].right].right].right) || boundary_indices.contains(g.r[g.r[g.r[i].right].right].right) );
                }else {
                    will_3right = false;
                }

                if (will_2down) {
                    will_3down = (g.r[g.r[g.r[i].down].down].down != -1) && ( !inside_indices.contains(g.r[g.r[g.r[i].down].down].down) || boundary_indices.contains(g.r[g.r[g.r[i].down].down].down) );
                }else {
                    will_3down = false;
                }

                if (will_2up) {
                    will_3up = (g.r[g.r[g.r[i].up].up].up != -1) && ( !inside_indices.contains(g.r[g.r[g.r[i].up].up].up) || boundary_indices.contains(g.r[g.r[g.r[i].up].up].up) );
                }else {
                    will_3up = false;
                }

                if (will_2front) {
                    will_3front = (g.r[g.r[g.r[i].front].front].front != -1) && ( !inside_indices.contains(g.r[g.r[g.r[i].front].front].front) || boundary_indices.contains(g.r[g.r[g.r[i].front].front].front) );
                }else {
                    will_3front = false;
                }

                if (will_2back) {
                    will_3back = (g.r[g.r[g.r[i].back].back].back != -1) && ( !inside_indices.contains(g.r[g.r[g.r[i].back].back].back) || boundary_indices.contains(g.r[g.r[g.r[i].back].back].back) );
                }else {
                    will_3back = false;
                }




                //std::cerr << "\t3\n";
                const bool can_cx = will_2right && will_2left;
                const bool can_fx = will_3right;
                const bool can_bx = will_3left;

                const bool can_cy = will_2up && will_2down;
                const bool can_fy = will_3up;
                const bool can_by = will_3down;

                const bool can_cz = will_2front && will_2back;
                const bool can_fz = will_3back;
                const bool can_bz = will_3front;

                //const bool can_cx
                const bool can_x = can_cx || can_fx || can_bx;
                const bool can_y = can_cy || can_fy || can_by;
                const bool can_z = can_cz || can_fz || can_bz;

                const bool can_deriv = can_x && can_y && can_z;
                //std::cerr << "\t4\n";


                if (!can_deriv) {
                    bool did_something = false;
                    //std::cerr << "can't deriv\n";
                    //moving the points in a direction
                    if (inside_indices.contains(g.r[i].left) && !inside_indices.contains(g.r[i].right)) {
                        boundary_indices.insert( g.r[i].right );
                        inside_indices.insert( g.r[i].right );

                        norms.add_point(g.r[i].right, norms.normal(i));
                        t_inds.add_point(g.r[i].right, t_inds.get_index(i));
                        m_points.add_point(g.r[i].right, m_points.get_point(i));
                        v_points.add_point(g.r[i].right, v_points.get_vel(i));

                        did_something = true;

                    } else if (!inside_indices.contains(g.r[i].left) && inside_indices.contains(g.r[i].right)) {
                        boundary_indices.insert( g.r[i].left );
                        inside_indices.insert( g.r[i].left );

                        norms.add_point(g.r[i].left, norms.normal(i));
                        t_inds.add_point(g.r[i].left, t_inds.get_index(i));
                        m_points.add_point(g.r[i].left, m_points.get_point(i));
                        v_points.add_point(g.r[i].left, v_points.get_vel(i));

                        did_something = true;
                    }

                    //std::cerr << "\t9\n";
                    if (inside_indices.contains(g.r[i].up) && !inside_indices.contains(g.r[i].down)) {
                        //std::cerr<<"1\n";
                        boundary_indices.insert( g.r[i].down );
                        inside_indices.insert( g.r[i].down );

                        norms.add_point(g.r[i].down, norms.normal(i));
                        t_inds.add_point(g.r[i].down, t_inds.get_index(i));
                        m_points.add_point(g.r[i].down, m_points.get_point(i));
                        v_points.add_point(g.r[i].down, v_points.get_vel(i));

                        did_something = true;

                    } else if (!inside_indices.contains(g.r[i].up) && inside_indices.contains(g.r[i].down)) {
                        /*std::cerr<<"2\n";
                        std::cerr << g.r[i].up << "\n";
                        std::cerr << g.r[i].down << "\n";

                        std::cerr << norms.normal(i) << "\n";
                        std::cerr << t_inds.get_index(i) << "\n";
                        std::cerr << m_points.get_point(i) << "\n";
                        std::cerr << v_points.get_vel(i) << "\n";*/

                        boundary_indices.insert( g.r[i].up );
                        inside_indices.insert( g.r[i].up );

                        norms.add_point(g.r[i].up, norms.normal(i));
                        t_inds.add_point(g.r[i].up, t_inds.get_index(i));
                        m_points.add_point(g.r[i].up, m_points.get_point(i));
                        v_points.add_point(g.r[i].up, v_points.get_vel(i));

                        did_something = true;
                    }
                    //std::cerr << "\t10\n";
                    if (inside_indices.contains(g.r[i].front) && !inside_indices.contains(g.r[i].back)) {
                        boundary_indices.insert( g.r[i].back );
                        inside_indices.insert( g.r[i].back );

                        norms.add_point(g.r[i].back, norms.normal(i));
                        t_inds.add_point(g.r[i].back, t_inds.get_index(i));
                        m_points.add_point(g.r[i].back, m_points.get_point(i));
                        v_points.add_point(g.r[i].back, v_points.get_vel(i));

                        did_something = true;

                    } else if (!inside_indices.contains(g.r[i].front) && inside_indices.contains(g.r[i].back)) {
                        boundary_indices.insert( g.r[i].front );
                        inside_indices.insert( g.r[i].front );

                        norms.add_point(g.r[i].front, norms.normal(i));
                        t_inds.add_point(g.r[i].front, t_inds.get_index(i));
                        m_points.add_point(g.r[i].front, m_points.get_point(i));
                        v_points.add_point(g.r[i].front, v_points.get_vel(i));

                        did_something = true;
                    }





                    //std::cerr << "\t11\n";


                    if (!can_x && inside_indices.contains(g.r[i].left) && inside_indices.contains(g.r[i].right) ) {
                        //try to pick the best direction to add the point
                        if (!inside_indices.contains(g.r[g.r[i].left].left)) {    //best dir left
                            boundary_indices.insert( g.r[i].left );
                            inside_indices.insert( g.r[i].left );

                            norms.add_point(g.r[i].left, norms.normal(i));
                            t_inds.add_point(g.r[i].left, t_inds.get_index(i));
                            m_points.add_point(g.r[i].left, m_points.get_point(i));
                            v_points.add_point(g.r[i].left, v_points.get_vel(i));

                            /*const auto &ind3 = g.r[g.r[g.r[i].left].left].left;

                            boundary_indices.insert( ind3 );
                            inside_indices.insert( ind3 );

                            norms.add_point(ind3, norms.normal(i));
                            t_inds.add_point(ind3, t_inds.get_index(i));
                            m_points.add_point(ind3, m_points.get_point(i));
                            v_points.add_point(ind3, v_points.get_vel(i));*/

                        } else if (!inside_indices.contains(g.r[g.r[i].right].right)){   //best dir right
                            boundary_indices.insert( g.r[i].right );
                            inside_indices.insert( g.r[i].right );

                            norms.add_point(g.r[i].right, norms.normal(i));
                            t_inds.add_point(g.r[i].right, t_inds.get_index(i));
                            m_points.add_point(g.r[i].right, m_points.get_point(i));
                            v_points.add_point(g.r[i].right, v_points.get_vel(i));

                            /*const auto &ind3 = g.r[g.r[g.r[i].right].right].right;

                            boundary_indices.insert( ind3 );
                            inside_indices.insert( ind3 );

                            norms.add_point(ind3, norms.normal(i));
                            t_inds.add_point(ind3, t_inds.get_index(i));
                            m_points.add_point(ind3, m_points.get_point(i));
                            v_points.add_point(ind3, v_points.get_vel(i));*/
                        } else {    //can't find best direction so just add both directions
                            boundary_indices.insert( g.r[i].right );
                            inside_indices.insert( g.r[i].right );

                            norms.add_point(g.r[i].right, norms.normal(i));
                            t_inds.add_point(g.r[i].right, t_inds.get_index(i));
                            m_points.add_point(g.r[i].right, m_points.get_point(i));
                            v_points.add_point(g.r[i].right, v_points.get_vel(i));


                            const auto &ind2_1 = g.r[g.r[i].right].right;
                            const auto &ind2_2 = g.r[g.r[i].left].left;
                            
                            boundary_indices.insert( ind2_1 );
                            inside_indices.insert( ind2_1 );

                            norms.add_point(ind2_1, norms.normal(i));
                            t_inds.add_point(ind2_1, t_inds.get_index(i));
                            m_points.add_point(ind2_1, m_points.get_point(i));
                            v_points.add_point(ind2_1, v_points.get_vel(i));


                            boundary_indices.insert( ind2_2 );
                            inside_indices.insert( ind2_2 );

                            norms.add_point(ind2_2, norms.normal(i));
                            t_inds.add_point(ind2_2, t_inds.get_index(i));
                            m_points.add_point(ind2_2, m_points.get_point(i));
                            v_points.add_point(ind2_2, v_points.get_vel(i));



                            boundary_indices.insert( g.r[i].left );
                            inside_indices.insert( g.r[i].left );

                            norms.add_point(g.r[i].left, norms.normal(i));
                            t_inds.add_point(g.r[i].left, t_inds.get_index(i));
                            m_points.add_point(g.r[i].left, m_points.get_point(i));
                            v_points.add_point(g.r[i].left, v_points.get_vel(i));
                        }
                    }
                    //std::cerr << "\t12\n";

                    if (!can_y && inside_indices.contains(g.r[i].down) && inside_indices.contains(g.r[i].up) ) {
                        //try to pick the best direction to add the point
                        if (!inside_indices.contains(g.r[g.r[i].down].down)) {    //best dir down

                            boundary_indices.insert( g.r[i].down );
                            inside_indices.insert( g.r[i].down );

                            norms.add_point(g.r[i].down, norms.normal(i));
                            t_inds.add_point(g.r[i].down, t_inds.get_index(i));
                            m_points.add_point(g.r[i].down, m_points.get_point(i));
                            v_points.add_point(g.r[i].down, v_points.get_vel(i));

                            /*boundary_indices.insert( g.r[g.r[i].down].down );
                            inside_indices.insert( g.r[g.r[i].down].down );

                            norms.add_point(g.r[g.r[i].down].down, norms.normal(i));
                            t_inds.add_point(g.r[g.r[i].down].down, t_inds.get_index(i));
                            m_points.add_point(g.r[g.r[i].down].down, m_points.get_point(i));
                            v_points.add_point(g.r[g.r[i].down].down, v_points.get_vel(i));*/

                            /*const auto &ind3 = g.r[g.r[g.r[i].down].down].down;

                            boundary_indices.insert( ind3 );
                            inside_indices.insert( ind3 );

                            norms.add_point(ind3, norms.normal(i));
                            t_inds.add_point(ind3, t_inds.get_index(i));
                            m_points.add_point(ind3, m_points.get_point(i));
                            v_points.add_point(ind3, v_points.get_vel(i));*/

                        } else if (!inside_indices.contains(g.r[g.r[i].up].up)){   //best dir up

                            boundary_indices.insert( g.r[i].up );
                            inside_indices.insert( g.r[i].up );

                            norms.add_point(g.r[i].up, norms.normal(i));
                            t_inds.add_point(g.r[i].up, t_inds.get_index(i));
                            m_points.add_point(g.r[i].up, m_points.get_point(i));
                            v_points.add_point(g.r[i].up, v_points.get_vel(i));

                            /*boundary_indices.insert( g.r[g.r[i].up].up );
                            inside_indices.insert( g.r[g.r[i].up].up );

                            norms.add_point(g.r[g.r[i].up].up, norms.normal(i));
                            t_inds.add_point(g.r[g.r[i].up].up, t_inds.get_index(i));
                            m_points.add_point(g.r[g.r[i].up].up, m_points.get_point(i));
                            v_points.add_point(g.r[g.r[i].up].up, v_points.get_vel(i));*/

                            /*const auto &ind3 = g.r[g.r[g.r[i].up].up].up;

                            boundary_indices.insert( ind3 );
                            inside_indices.insert( ind3 );

                            norms.add_point(ind3, norms.normal(i));
                            t_inds.add_point(ind3, t_inds.get_index(i));
                            m_points.add_point(ind3, m_points.get_point(i));
                            v_points.add_point(ind3, v_points.get_vel(i));*/
                        } else {    //can't find best direction so just add both directions

                            boundary_indices.insert( g.r[i].up );
                            inside_indices.insert( g.r[i].up );

                            norms.add_point(g.r[i].up, norms.normal(i));
                            t_inds.add_point(g.r[i].up, t_inds.get_index(i));
                            m_points.add_point(g.r[i].up, m_points.get_point(i));
                            v_points.add_point(g.r[i].up, v_points.get_vel(i));


                            boundary_indices.insert( g.r[i].down );
                            inside_indices.insert( g.r[i].down );

                            norms.add_point(g.r[i].down, norms.normal(i));
                            t_inds.add_point(g.r[i].down, t_inds.get_index(i));
                            m_points.add_point(g.r[i].down, m_points.get_point(i));
                            v_points.add_point(g.r[i].down, v_points.get_vel(i));


                            const auto &ind2_1 = g.r[g.r[i].up].up;
                            const auto &ind2_2 = g.r[g.r[i].down].down;

                            boundary_indices.insert( ind2_1 );
                            inside_indices.insert( ind2_1 );

                            norms.add_point(ind2_1, norms.normal(i));
                            t_inds.add_point(ind2_1, t_inds.get_index(i));
                            m_points.add_point(ind2_1, m_points.get_point(i));
                            v_points.add_point(ind2_1, v_points.get_vel(i));


                            boundary_indices.insert( ind2_2 );
                            inside_indices.insert( ind2_2 );

                            norms.add_point(ind2_2, norms.normal(i));
                            t_inds.add_point(ind2_2, t_inds.get_index(i));
                            m_points.add_point(ind2_2, m_points.get_point(i));
                            v_points.add_point(ind2_2, v_points.get_vel(i));

                        }
                    }


                    //std::cerr << "\t13\n";
                    if (!can_z && inside_indices.contains(g.r[i].back) && inside_indices.contains(g.r[i].front) ) {
                        //try to pick the best direction to add the point
                        if (!inside_indices.contains(g.r[g.r[i].back].back)) {    //best dir back
                            //std::cerr<<"1\n";
                            boundary_indices.insert( g.r[i].back );
                            inside_indices.insert( g.r[i].back );

                            norms.add_point(g.r[i].back, norms.normal(i));
                            t_inds.add_point(g.r[i].back, t_inds.get_index(i));
                            m_points.add_point(g.r[i].back, m_points.get_point(i));
                            v_points.add_point(g.r[i].back, v_points.get_vel(i));

                            /*boundary_indices.insert( g.r[g.r[i].back].back );
                            inside_indices.insert( g.r[g.r[i].back].back );

                            norms.add_point(g.r[g.r[i].back].back, norms.normal(i));
                            t_inds.add_point(g.r[g.r[i].back].back, t_inds.get_index(i));
                            m_points.add_point(g.r[g.r[i].back].back, m_points.get_point(i));
                            v_points.add_point(g.r[g.r[i].back].back, v_points.get_vel(i));*/

                            /*const auto &ind3 = g.r[g.r[g.r[i].back].back].back;

                            boundary_indices.insert( ind3 );
                            inside_indices.insert( ind3 );

                            norms.add_point(ind3, norms.normal(i));
                            t_inds.add_point(ind3, t_inds.get_index(i));
                            m_points.add_point(ind3, m_points.get_point(i));
                            v_points.add_point(ind3, v_points.get_vel(i));*/

                        } else if (!inside_indices.contains(g.r[g.r[i].front].front)){   //best dir front
                            //std::cerr<<"2\n";
                            boundary_indices.insert( g.r[i].front );
                            inside_indices.insert( g.r[i].front );

                            norms.add_point(g.r[i].front, norms.normal(i));
                            t_inds.add_point(g.r[i].front, t_inds.get_index(i));
                            m_points.add_point(g.r[i].front, m_points.get_point(i));
                            v_points.add_point(g.r[i].front, v_points.get_vel(i));

                            /*boundary_indices.insert( g.r[g.r[i].front].front );
                            inside_indices.insert( g.r[g.r[i].front].front );

                            norms.add_point(g.r[g.r[i].front].front, norms.normal(i));
                            t_inds.add_point(g.r[g.r[i].front].front, t_inds.get_index(i));
                            m_points.add_point(g.r[g.r[i].front].front, m_points.get_point(i));
                            v_points.add_point(g.r[g.r[i].front].front, v_points.get_vel(i));*/

                            /*const auto &ind3 = g.r[g.r[g.r[i].front].front].front;

                            boundary_indices.insert( ind3 );
                            inside_indices.insert( ind3 );

                            norms.add_point(ind3, norms.normal(i));
                            t_inds.add_point(ind3, t_inds.get_index(i));
                            m_points.add_point(ind3, m_points.get_point(i));
                            v_points.add_point(ind3, v_points.get_vel(i));*/
                        } else {    //can't find best direction so just add both directions
                            /*std::cerr<<"3\n";
                            std::cerr << g.r[i].front << "\n";
                            std::cerr << g.r[i].back << "\n";

                            std::cerr << norms.normal(i) << "\n";
                            std::cerr << t_inds.get_index(i) << "\n";
                            std::cerr << m_points.get_point(i) << "\n";
                            std::cerr << v_points.get_vel(i) << "\n";*/

                            boundary_indices.insert( g.r[i].front );
                            inside_indices.insert( g.r[i].front );

                            norms.add_point(g.r[i].front, norms.normal(i));
                            t_inds.add_point(g.r[i].front, t_inds.get_index(i));
                            m_points.add_point(g.r[i].front, m_points.get_point(i));
                            v_points.add_point(g.r[i].front, v_points.get_vel(i));


                            boundary_indices.insert( g.r[i].back );
                            inside_indices.insert( g.r[i].back );

                            norms.add_point(g.r[i].back, norms.normal(i));
                            t_inds.add_point(g.r[i].back, t_inds.get_index(i));
                            m_points.add_point(g.r[i].back, m_points.get_point(i));
                            v_points.add_point(g.r[i].back, v_points.get_vel(i));


                            const auto &ind2_1 = g.r[g.r[i].back].back;
                            const auto &ind2_2 = g.r[g.r[i].front].front;

                            boundary_indices.insert( ind2_1 );
                            inside_indices.insert( ind2_1 );

                            norms.add_point(ind2_1, norms.normal(i));
                            t_inds.add_point(ind2_1, t_inds.get_index(i));
                            m_points.add_point(ind2_1, m_points.get_point(i));
                            v_points.add_point(ind2_1, v_points.get_vel(i));


                            boundary_indices.insert( ind2_2 );
                            inside_indices.insert( ind2_2 );

                            norms.add_point(ind2_2, norms.normal(i));
                            t_inds.add_point(ind2_2, t_inds.get_index(i));
                            m_points.add_point(ind2_2, m_points.get_point(i));
                            v_points.add_point(ind2_2, v_points.get_vel(i));
                        }
                    }
                   // std::cerr<<"4\n";


                    if (did_something) {
                        boundary_indices.erase(i);
                        norms.erase(i);
                        t_inds.erase(i);
                        m_points.erase(i);
                        v_points.erase(i);
                    }

                    //std::cerr << "\t14\n";
                    inside_indices.insert(i);
                    /*boundary_indices.erase(i);
                    norms.erase(i);
                    t_inds.erase(i);
                    m_points.erase(i);
                    v_points.erase(i);*/

                    was_update = true;
                    //std::cerr << "\t15\n";
                    //std::cerr << "point " << i << " needed updating\n";
                    //std::cerr << "x: " << can_x << " y: " << can_y << " z: " << can_z << "\n";

                }

            }
        }
        //std::cerr << "==========================================================\n";
    }

    //don't want boundary points infron of other boundary points
    /*std::unordered_set<unsigned> indices_to_erase;
    for (unsigned i = 0; i < g.size(); i++) {
        if (boundary_indices.contains(i)) {
            //get direction where normal is biggest
            const auto &n = norms.normal(i);


            constexpr double mc = 0.8;  //minimum a normal vector needs to be in a direction for that direction to be considered

            if (std::abs(n.x()) > mc) {
                if (boundary_indices.contains( g.get_move_ind(i, 1,0,0) )  &&  boundary_indices.contains( g.get_move_ind(i, -1,0,0) )) {
                    indices_to_erase.insert(i);
                }
            }
            if (std::abs(n.y()) > mc) {
                if (boundary_indices.contains( g.get_move_ind(i, 0,1,0) )  &&  boundary_indices.contains( g.get_move_ind(i, 0,-1,0) )) {
                    indices_to_erase.insert(i);
                }
            }
            if (std::abs(n.z()) > mc) {
                if (boundary_indices.contains( g.get_move_ind(i, 0,0,1) )  &&  boundary_indices.contains( g.get_move_ind(i, 0,0,-1) )) {
                    indices_to_erase.insert(i);
                }

            }

        }
    }

    for (const auto& err : indices_to_erase) {
        boundary_indices.erase(err);
        norms.erase(err);
        t_inds.erase(err);
        m_points.erase(err);
        v_points.erase(err);
    }*/

#ifdef REMOVE_INSIDE_DLOG
    std::cerr << "removing points\n";
#endif

    unsigned index_counter = 0;

    for (unsigned i = 0; i < g.x.size(); i++) {

        if (inside_indices.contains(i) && !boundary_indices.contains(i) || g.r[i].is_edge() ) { // is_edge is for removing edge points at walls
            old_new.insert({i, -1});
        } else {
            x.push_back(g.x[i]);
            y.push_back(g.y[i]);
            z.push_back(g.z[i]);
            r.push_back(g.r[i]);

            if (boundary_indices.contains(i)) {
                if (inside_indices.contains(r[index_counter].left) && !boundary_indices.contains(r[index_counter].left)) { r[index_counter].left = -1; }
                if (inside_indices.contains(r[index_counter].right) && !boundary_indices.contains(r[index_counter].right)) { r[index_counter].right = -1; }
                if (inside_indices.contains(r[index_counter].down) && !boundary_indices.contains(r[index_counter].down)) { r[index_counter].down = -1; }
                if (inside_indices.contains(r[index_counter].up) && !boundary_indices.contains(r[index_counter].up)) { r[index_counter].up = -1; }
                if (inside_indices.contains(r[index_counter].front) && !boundary_indices.contains(r[index_counter].front)) { r[index_counter].front = -1; }
                if (inside_indices.contains(r[index_counter].back) && !boundary_indices.contains(r[index_counter].back)) { r[index_counter].back = -1; }

                /*const auto &n = norms.normal(i);
                //getting looking directions
                // - look in the direction where the normal is bigger than 1/sqrt(3)
                // - check to see if in the reverse direction there is a boundary
                // - if there is a boundary point, treat that point as off the grid
                if (std::abs(n.x()) > 1/sqrt(3)) {
                    if (n.x() > 0) {
                        if (boundary_indices.contains(r[index_counter].left)) { r[index_counter].left = -1; }
                    } else {
                        if (boundary_indices.contains(r[index_counter].right)) { r[index_counter].right = -1; }
                    }
                }
                if (std::abs(n.y()) > 1/sqrt(3)) {
                    if (n.y() > 0) {
                        if (boundary_indices.contains(r[index_counter].down)) { r[index_counter].down = -1; }
                    } else {
                        if (boundary_indices.contains(r[index_counter].up)) { r[index_counter].up = -1; }
                    }
                }
                if (std::abs(n.z()) > 1/sqrt(3)) {
                    if (n.z() > 0) {
                        if (boundary_indices.contains(r[index_counter].front)) { r[index_counter].front = -1; }
                    } else {
                        if (boundary_indices.contains(r[index_counter].back)) { r[index_counter].back = -1; }
                    }
                }*/



                /*if (boundary_indices.contains(r[index_counter].left)) { r[index_counter].left = -1; }
                if (boundary_indices.contains(r[index_counter].right)) { r[index_counter].right = -1; }
                if (boundary_indices.contains(r[index_counter].down)) { r[index_counter].down = -1; }
                if (boundary_indices.contains(r[index_counter].up)) { r[index_counter].up = -1; }
                if (boundary_indices.contains(r[index_counter].front)) { r[index_counter].front = -1; }
                if (boundary_indices.contains(r[index_counter].back)) { r[index_counter].back = -1; }


                if (inside_indices.contains(r[index_counter].left) ) { r[index_counter].left = -1; }
                if (inside_indices.contains(r[index_counter].right) ) { r[index_counter].right = -1; }
                if (inside_indices.contains(r[index_counter].down) ) { r[index_counter].down = -1; }
                if (inside_indices.contains(r[index_counter].up) ) { r[index_counter].up = -1; }
                if (inside_indices.contains(r[index_counter].front) ) { r[index_counter].front = -1; }
                if (inside_indices.contains(r[index_counter].back) ) { r[index_counter].back = -1; }*/
                //update_r_point(r[index_counter], boundary_indices, inside_indices);
            }

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

    g.set_plotting_points();


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

            #ifdef REMOVE_INSIDE_DLOG
    std::cerr << "recalibrating boundary points\n";
            #endif
    for (const auto & i : boundary_indices) {
        g.boundary_indices.insert( old_new[i]  );
    }

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

            #ifdef REMOVE_INSIDE_DLOG
    std::cerr << "done\n";
            #endif

            #ifdef CHECK_GRID_RESULTS
            #ifdef REMOVE_INSIDE_DLOG
    std::cerr << "checking results\n";
            #endif
    //checking that all boundary points have a normal
    for (unsigned i = 0; i < g.r.size(); i++) {
        if (g.is_boundary(i) && g.off_walls(i) && !norms.contains(i)) {
            std::cerr << "boundary point at index " << i << " does not have a normal\n";
            std::cerr << "\tThis is at (" << g.x[i] << " " << g.y[i] << " " << g.z[i] << ")\n";
        }
    }
    //checking that all boundary points have a normal
    for (unsigned i = 0; i < g.r.size(); i++) {
        if (g.is_boundary(i) && g.off_walls(i) && !v_points.contains(i)) {
            std::cerr << "boundary point at index " << i << " does not have a velocity\n";
            std::cerr << "\tThis is at (" << g.x[i] << " " << g.y[i] << " " << g.z[i] << ")\n";
        }
    }
    //checking that all boundary points are off wall
    for (const auto &v : v_points.m) {
        if (!g.off_walls(v.first)) {
            const auto i = v.first;
            std::cerr << "boundary point " << i << " is not off walls\n";
            std::cerr << "\tThis is at (" << g.x[i] << " " << g.y[i] << " " << g.z[i] << ")\n";

            const auto dist_left = dist_to_plane( {g.x[i], g.y[i], g.z[i]}, g.edge1, g.edge5, g.edge7 );
            const auto dist_right = dist_to_plane( {g.x[i], g.y[i], g.z[i]}, g.edge2, g.edge4, g.edge8 );
            const auto dist_down = dist_to_plane( {g.x[i], g.y[i], g.z[i]}, g.edge1, g.edge2, g.edge5 );
            const auto dist_up = dist_to_plane( {g.x[i], g.y[i], g.z[i]}, g.edge4, g.edge8, g.edge7 );
            const auto dist_front = dist_to_plane( {g.x[i], g.y[i], g.z[i]}, g.edge1, g.edge2, g.edge4 );
            const auto dist_back = dist_to_plane( {g.x[i], g.y[i], g.z[i]}, g.edge5, g.edge7, g.edge8 );

            std::cerr << "\tl:" << dist_left << " r:" << dist_right << " d:" << dist_down << " u:" << dist_up << " f:" << dist_front << " b:" << dist_back << "\n";
            std::cerr << "\tdx: "<< g.dx << " dy: " << g.dy << " dz:" << g.dz << "\n";

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
