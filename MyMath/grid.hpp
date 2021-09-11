//
// Created by jacob on 11/9/21.
//

#ifndef CODE_GRID_HPP
#define CODE_GRID_HPP

#include <vector>

//stores the neighbours of a grid point
// -1 if no neighbour
struct grid_relation {
    int left{}, right{}, down{}, up{}, front{}, back{};
};

//grids must be axis aligned
struct grid {
    std::vector<double> x, y, z;
    double dx, dy, dz;
    std::vector<grid_relation> r{};


    grid() : dx{}, dy{}, dz{} {};
    grid(const std::vector<double> &x_, const std::vector<double> &y_, const std::vector<double> &z_, const double dx_, const double dy_, const double dz_)
        : x(x_), y(y_), z(z_), dx(dx_), dy(dy_), dz(dz_) {
#ifndef NDEBUG
        bool err = false;
        if (x_.size() != y_.size()) {
            std::cerr << "x and y need to be the same size\n";
            err = true;
        }
        if (y_.size() != z_.size()) {
            std::cerr <<"y and z need to be the same size\n";
            err = true;
        }
        if (x_.size() != z_.size()) {
            std::cerr << "x and z need to be the same size\n";
            err = true;
        }
        if (err) {
            std::cerr << "\tx.size = " << x_.size() << "\ty.size = " << y_.size() << "\tz.size = " << z_.size() << "\n";
        }
#endif
    }
};

#endif //CODE_GRID_HPP
