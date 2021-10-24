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

//#define BOUNDARY_CONDITIONS_DLOG


struct boundary_conditions {
    size_t no_inside_mesh = 0;
    grid global_grid;
    boundary_normals_grouped norms; //need to differentiate between wall norms and mesh norms. Still want to keep norms to not have to change a lot of code.
    boundary_normals mesh_norms;
    wall_normals wall_norms;
    mesh_points m_points;
    vel_points v_points;
    tri_inds t_inds;
    std::map<unsigned, int> old_new;    //conversion between the old indices and the new (points removed) indices

    triangle_mesh tm;


    boundary_conditions() = delete;
    boundary_conditions(const mesh *m, const double dx, const double dy, const double dz, const double Wx, const double Wy, const double Wz, const unsigned sx, const unsigned sy, const unsigned sz, const double minx = 0, const double miny = 0, const double minz = 0)
        : tm(m), norms(mesh_norms, wall_norms) {
#ifdef BOUNDARY_CONDITIONS_DLOG
        std::cerr << "making entire grid\n";
#endif
        make_entire_grid(global_grid, Wx, Wy, Wz, dx, dy, dz, sx, sy, sz, minx, miny, minz);


#ifdef BOUNDARY_CONDITIONS_DLOG
        std::cerr << "removing inside points\n";
#endif
        remove_inside_boundary_unif(global_grid, tm, *m, mesh_norms, m_points, v_points, old_new, t_inds);

#ifdef BOUNDARY_CONDITIONS_DLOG
        std::cerr << "setting boundary normals\n";
#endif
        wall_norms = wall_normals(no_wall_points());
#ifdef BOUNDARY_CONDITIONS_DLOG
        std::cerr << "creating wall normals\n";
#endif
        create_wall_normals();
    }

    //should never be used in real flow, only used for testing derivatives
    boundary_conditions(const double dx, const double dy, const double dz, const double Wx, const double Wy, const double Wz, const unsigned sx, const unsigned sy, const unsigned sz, const double minx = 0, const double miny = 0, const double minz = 0)
        : tm{}, norms(mesh_norms, wall_norms) {
        make_entire_grid(global_grid, Wx, Wy, Wz, dx, dy, dz, sx, sy, sz, minx, miny, minz);

        wall_norms = wall_normals(no_wall_points());
        create_wall_normals();

        global_grid.set_plotting_points();
    }

    void update(const vec3& w, const vec3& v, const vec3 &c_o_m, const double dt) {
        tm.update();
        m_points.update(w, v, c_o_m, dt);
        v_points.update(t_inds, m_points);
    }

    [[nodiscard]] auto size() const {
        return m_points.size();
    }

    unsigned no_wall_points() const {
        const auto dims = global_grid.no_points_unif;
        const auto N = static_cast<unsigned>(dims.x()-1);
        const auto M = static_cast<unsigned>(dims.y()-1);
        const auto P = static_cast<unsigned>(dims.z()-1);
        return 2*(N+1)*(P+1) + 2*(N-1)*(M+1) + 2*(M-1)*(P-1);
    }

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

    void DEBUG_write_normal_vectors_at_x(const double xp) const {
        std::ofstream output("../DEBUG/normal_vectors_x.txt");
        if (output.is_open()) {
            const auto inds = global_grid.get_some_x_inds(xp);
            for (const auto ind : inds) {
                output << global_grid[ind].y() << " " << global_grid[ind].z() << " ";
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



    //void update_velocity_wall_BC();

private:
    void create_wall_normals();
};

void boundary_conditions::create_wall_normals() {
    const auto dims = global_grid.no_points_unif;
    const auto N = static_cast<unsigned>(dims.x()-1);
    const auto M = static_cast<unsigned>(dims.y()-1);
    const auto P = static_cast<unsigned>(dims.z()-1);

    for (unsigned i = 1; i < N; i++) {
        for (unsigned k = 1; k < P; k++) {
            const auto ind1 = old_new[global_grid.convert_indices_unif(vec3(i,0,k) )];
            wall_norms.add_point(ind1, vec3(0,1,0));
            const auto ind2 = old_new[global_grid.convert_indices_unif(vec3(i,M,k) ) ];
            wall_norms.add_point(ind2, vec3(0,-1,0));
        }
    }
    for (unsigned i = 1; i < N; i++) {
        for (unsigned j = 1; j < M; j++) {
            const auto ind1 = old_new[global_grid.convert_indices_unif(vec3(i,j,0) )];
            wall_norms.add_point(ind1, vec3(0,0,1));
            const auto ind2 = old_new[global_grid.convert_indices_unif(vec3(i,j,P) )];
            wall_norms.add_point(ind2, vec3(0,0,-1));
        }
    }
    for (unsigned j = 1; j < M; j++) {
        for (unsigned k = 1; k < P; k++) {
            const auto ind1 = old_new[global_grid.convert_indices_unif(vec3(0,j,k) )];
            wall_norms.add_point(ind1,vec3(1,0,0));
            const auto ind2 = old_new[global_grid.convert_indices_unif(vec3(N,j,k) )];
            wall_norms.add_point(ind2,vec3(-1,0,0));
        }
    }


}


#endif //CODE_BOUNDARY_CONDITIONS_HPP
