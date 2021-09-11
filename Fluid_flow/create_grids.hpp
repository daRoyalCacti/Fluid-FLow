//
// Created by jacob on 11/9/21.
//

#ifndef CODE_CREATE_GRIDS_HPP
#define CODE_CREATE_GRIDS_HPP

#include "../MyMath/grid.hpp"
#include <vector>

void make_enitre_gid(grid &g, const double Wx, const double Wy, const double Wz, const double dx, const double dy, const double dz) {
    g.x.resize(static_cast<unsigned long>(Wx/dx));
    g.y.resize(static_cast<unsigned long>(Wy/dy));
    g.z.resize(static_cast<unsigned long>(Wz/dz));

    unsigned counter = 0;
    for (unsigned i = 0; i < g.x.size(); i++) {
        for (unsigned j = 0; j < g.y.size(); j++) {
            for (unsigned k = 0; k < g.z.size(); k++) {
                g.x[counter] = i*dx;
                g.y[counter] = j*dy;
                g.z[counter] = k*dz;

                counter++;
            }
        }
    }

}

#endif //CODE_CREATE_GRIDS_HPP
