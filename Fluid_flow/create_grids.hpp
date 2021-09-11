//
// Created by jacob on 11/9/21.
//

#ifndef CODE_CREATE_GRIDS_HPP
#define CODE_CREATE_GRIDS_HPP

#include "../MyMath/grid.hpp"
#include <vector>

void make_enitre_gid(grid &g, const double Wx, const double Wy, const double Wz, const double dx, const double dy, const double dz) {
    const auto sx = static_cast<unsigned long>(Wx/dx);
    const auto sy = static_cast<unsigned long>(Wy/dy);
    const auto sz = static_cast<unsigned long>(Wz/dz);

    const auto s = static_cast<unsigned long>(sx*sy*sz);
    g.x.resize(s);
    g.y.resize(s);
    g.z.resize(s);

    unsigned counter = 0;
    for (unsigned i = 0; i < sx; i++) {
        for (unsigned j = 0; j < sy; j++) {
            for (unsigned k = 0; k < sz; k++) {
                g.x[counter] = i*dx;
                g.y[counter] = j*dy;
                g.z[counter] = k*dz;

                counter++;
            }
        }
    }

}

#endif //CODE_CREATE_GRIDS_HPP
