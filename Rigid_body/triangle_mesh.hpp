//
// Created by jacob on 28/8/21.
//

#ifndef CODE_TRIANGLE_MESH_HPP
#define CODE_TRIANGLE_MESH_HPP

#include "triangle.hpp"

struct triangle_mesh {
    std::vector<triangle> tris;

    triangle_mesh() = delete;
    explicit triangle_mesh(const mesh *m) {
        tris.resize(m->indices.size() / 3);
        for (unsigned i = 0; i < tris.size(); i++) {
            tris[i] = triangle( &m->get_vertice_index(3*i), &m->get_vertice_index(3*i+1), &m->get_vertice_index(3*i+2),
                                &m->get_normal_index(3*i), &m->get_normal_index(3*i+1), &m->get_normal_index(3*i+2));
        }

    }

};

#endif //CODE_TRIANGLE_MESH_HPP
