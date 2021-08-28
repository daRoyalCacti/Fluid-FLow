//
// Created by jacob on 28/8/21.
//

#ifndef CODE_TRIANGLE_MESH_HPP
#define CODE_TRIANGLE_MESH_HPP

#include "ray.hpp"
#include "triangle.hpp"

template <class T>
class no_init_alloc
    : public std::allocator<T>
    {
    public:
        using std::allocator<T>::allocator;

        template <class U, class... Args> void construct(U*, Args&&...) {}
    };

struct triangle_mesh {
    std::vector<triangle, no_init_alloc<triangle>> tris;    //custom allocator because vector is filled immediately after resizing it

    triangle_mesh() = delete;
    explicit triangle_mesh(const mesh *m) {
        tris.resize(m->indices.size() / 3);
        for (unsigned i = 0; i < tris.size(); i++) {
            tris[i] = triangle( &m->get_vertice_index(3*i), &m->get_vertice_index(3*i+1), &m->get_vertice_index(3*i+2),
                                &m->get_normal_index(3*i), &m->get_normal_index(3*i+1), &m->get_normal_index(3*i+2));
        }

    }

    //returns if ray collided with mesh
    bool get_collision_points(const ray &r, vec3 &col1, vec3 &col2) {
        bool has_hit = false;
        bool hit_twice = false;
        double time;
        for (const auto & t : tris) {
            if (t.hit_time(r, time)) {
                if (has_hit) {
#ifndef NDEBUG
                    if (hit_twice) {
                        std::cerr << "Ray hit triangle mesh more than twice --- this is not implemented\n";
                    }
#endif
                    hit_twice = true;
                    col2 = r.at(time);
                } else {
                    has_hit = true;
                    col1 = r.at(time);
                }
            }
        }

        //if only hit once
        // - should rarely happen - only when ray just grazes past a triangle
        if (has_hit && !hit_twice) {
            col2 = col1;    //say both collisions happened at the same place
        }

        return has_hit;
    }


};

#endif //CODE_TRIANGLE_MESH_HPP
