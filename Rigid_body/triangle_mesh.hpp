//
// Created by jacob on 28/8/21.
//

#ifndef CODE_TRIANGLE_MESH_HPP
#define CODE_TRIANGLE_MESH_HPP

#include "ray.hpp"
#include "triangle.hpp"
#include "mesh.hpp"
#include <map>
#include <algorithm>


struct triangle_mesh {
    template <class T>
    class no_init_alloc
            : public std::allocator<T>
            {
            public:
                using std::allocator<T>::allocator;

                template <class U, class... Args> void construct(U*, Args&&...) {}
            };

    std::vector<triangle, no_init_alloc<triangle>> tris;    //custom allocator because vector is filled immediately after resizing it

    triangle_mesh() = default;
    explicit triangle_mesh(const mesh *m) {
        /*tris.resize(m->indices.size() / 3);
        for (unsigned i = 0; i < tris.size(); i++) {
            tris[i] = triangle( &m->get_vertice_index(3*i), &m->get_vertice_index(3*i+1), &m->get_vertice_index(3*i+2),
                                &m->get_normal_index(3*i), &m->get_normal_index(3*i+1), &m->get_normal_index(3*i+2),
                                &m->get_velocity_index(3*i), &m->get_velocity_index(3*i+1), &m->get_velocity_index(3*i+2),
                                &(m->v), &(m->w), &(m->pos_cm));
        }*/
        tris.resize(m->indices.size() / 3);
        for (unsigned i = 0; i < m->indices.size(); i+=3) {
            tris[i/3] = triangle( &m->get_vertice_index(i), &m->get_vertice_index(i+1), &m->get_vertice_index(i+2),
                                &m->get_normal_index(i), &m->get_normal_index(i+1), &m->get_normal_index(i+2),
                                &m->get_velocity_index(i), &m->get_velocity_index(i+1), &m->get_velocity_index(i+2),
                                &(m->v), &(m->w), &(m->pos_cm));
        }

    }

    /*struct d_vec3 {
        vec3 v1, v2, v3;
        const triangle* v4;
    };*/

    //returns if ray collided with mesh
    //[[nodiscard]] std::map<double, d_vec3> get_collision_points(const ray &r) const {
    [[nodiscard]] std::map<double, const triangle*> get_collision_points(const ray &r) const {
        //this does not work if the ray hits the mesh an odd number of times
        // - should only happen when a ray just barely grazes past the mesh
        // - note that hits is stored in order of hit time
        std::map<double, const triangle*> hits; //stores hit time and hit pos
        double time;
        for (const triangle & t : tris) {
            if (t.hit_time(r, time)) {
                //time, collision point, normal at collision point, velocity at collision point
                //hits.insert({ time, d_vec3{r.at(time), t.get_normal(r.at(time)), t.get_velocity(r.at(time)), &t}  });
                hits.insert({ time, &t  });
            }
        }

#ifndef NDEBUG
        if ( (hits.size() % 2) != 0 ) {
            std::cerr << "ray hit triangle mesh an odd number of times\n";
        }
#endif

        return hits;
    }

    void update() {
        for (auto &t : tris) {
            t.update();
        }
    }


};

#endif //CODE_TRIANGLE_MESH_HPP
