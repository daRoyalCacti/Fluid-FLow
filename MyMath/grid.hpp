//
// Created by jacob on 11/9/21.
//

#ifndef CODE_GRID_HPP
#define CODE_GRID_HPP

#include <vector>
#include <limits>

//stores the neighbours of a grid point
// -1 if no neighbour
struct grid_relation {
    int left{}, right{}, down{}, up{}, front{}, back{};
};

//grids must be axis aligned
// grid[i] corresponds to the interval x[i]+dx, y[i]+dy, z[i]+dz
struct grid {
    std::vector<double> x, y, z;
    double dx, dy, dz;
    std::vector<grid_relation> r{};
    vec3 mins, maxs;    //min = (minx, miny, minz)


    grid() : dx{}, dy{}, dz{} {};
    grid(const std::vector<double> &x_, const std::vector<double> &y_, const std::vector<double> &z_, const double dx_, const double dy_, const double dz_,
         const double minx, const double miny, const double minz, const double maxx, const double maxy, const double maxz)
        : x(x_), y(y_), z(z_), dx(dx_), dy(dy_), dz(dz_), mins(minx, miny, minz), maxs(maxx, maxy, maxz) {
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

    [[nodiscard]] constexpr inline auto get_ind(const vec3& p) const noexcept {
        //finds the index of the grid that p is in
        constexpr double inf = std::numeric_limits<double>::max();
        double min_x = inf, min_y = inf, min_z = inf;

        constexpr unsigned long init_ind = std::numeric_limits<unsigned long>::max(); //some stupid number
        unsigned long best_ind = init_ind;

        for (unsigned long i = 0; i < x.size(); i++) {
            const double dist_x = p.x() - x[i];
            if (dist_x < 0) {   //because the grid ranges from x[i],x[i]+dx being to the left of x[i] means the point is not in the grid
                continue;
            }
            const double dist_y = p.y() - y[i];
            if (dist_y < 0) {
                continue;
            }
            const double dist_z = p.z() - z[i];
            if (dist_z < 0) {
                continue;
            }

            if (dist_x < min_x || dist_y < min_y || dist_z < min_z) {
                min_x = dist_x;
                min_y = dist_y;
                min_z = dist_z;
                best_ind = i;
            }

        }

#ifndef NDEBUG
        if (best_ind == init_ind) {
            std::cerr << "point is to the left of the entire grid\n\tpoint = " << p << "\n";
        } else if (p.x() > x[best_ind]+dx || p.y() > y[best_ind]+dy || p.z() > z[best_ind] + dz ) { //must be an else if because x[best_ind] will error otherwise
            std::cerr << "point is outside of grid\n\tpoint = " << p << "\n";
        }
#endif

        return best_ind;

    }

};

#endif //CODE_GRID_HPP
