//
// Created by jacob on 11/9/21.
//

#ifndef CODE_GRID_HPP
#define CODE_GRID_HPP

#include <vector>
#include <limits>
#include <set>
#include "../MyMath/vec3.hpp"

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
    vec3 no_points_unif;    //the total number of points across each axis assuming the grid is boring


    vec3 operator[](unsigned i) const noexcept { return { x[i], y[i], z[i] }; }

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

    void update_grid(double x_off, double y_off, double z_off) noexcept {
        for (auto & x_ : x) {
            x_ += x_off;
        }
        for (auto & y_ : y) {
            y_ += y_off;
        }
        for (auto & z_ : z) {
            z_ += z_off;
        }
    }

    [[nodiscard]] auto size() const {
        return x.size();
    }

    [[nodiscard]] inline auto get_ind(const vec3& p) const noexcept {
        //finds the index of the grid that p is in
        for (unsigned long i = 0; i < x.size(); i++) {
            if (  x[i] <= p.x() && x[i]+dx >= p.x() //within the x bounds of the grid
            &&  y[i] <= p.y() && y[i]+dy >= p.y() //within the y bounds
            && z[i] <= p.z() && z[i]+dz >= p.z()) { //within the z bounds
                return i;
            }
        }

#ifndef NDEBUG
        std::cerr << "index was not found for the point\n";
#endif
        return std::numeric_limits<unsigned long>::max();   //some stupid value


    }

    [[nodiscard]] inline auto get_inds(vec3 p1, vec3 p2) const noexcept {
        //find the axis to look along
        unsigned axis1 = 4, axis2, axis3;
        const std::vector<double>* arr1, *arr2, *arr3;
        double d1, d2, d3;
        if (p1.x() != p2.x()) { //ray moving in x direction
            arr1 = &x;
            arr2 = &y;
            arr3 = &z;
            axis1 = 0;
            axis2 = 1;
            axis3 = 2;
            d1 = dx;
            d2 = dy;
            d3 = dz;
        } else if (p1.y() != p2.y()) {
            arr1 = &y;
            arr2 = &x;
            arr3 = &z;
            axis1 = 1;
            axis2 = 0;
            axis3 = 2;
            d1 = dy;
            d2 = dx;
            d3 = dz;
        } else if (p1.z() != p2.z()) {
            arr1 = &z;
            arr2 = &y;
            arr3 = &x;
            axis1 = 2;
            axis2 = 1;
            axis3 = 0;
            d1 = dz;
            d2 = dy;
            d3 = dx;
        }
#ifndef NDEBUG
        if (axis1 == 4) {
            std::cerr << "p1 cannot be equal to p2\n\tp1 = " << p1 << "\tp2 = " << p2 << "\n";
        }
        if (p1[axis2] != p2[axis2]) {
            std::cerr << "p1 and p2 cannot differ along 2 dimensions\n\tp1 = " << p1 << "\tp2 = " << p2 << "\n";
        }
        if (p1[axis3] != p2[axis3]) {
            std::cerr << "p1 and p2 cannot differ along 2 dimensions\n\tp1 = " << p1 << "\tp2 = " << p2 << "\n";
        }
#endif

        //ensuring p1 is to the left of p2
        if (p1[axis1] > p2[axis1]) {
            auto tmp = p1;
            p1 = p2;
            p2 = tmp;
        }



        std::vector<unsigned long> ret_vec;
        for (unsigned long i = 0; i < x.size(); i++) {
            if ( (*arr3)[i] < p1[axis3] && (*arr3)[i]+d3 > p2[axis3] ) {  //finding the points that lie on the same line as p1 and p2
                if ( (*arr2)[i] < p1[axis2] && (*arr2)[i]+d2 > p2[axis2] ) { //points that lie between p1 and p2 along axis 2
                    //the dimension the vectors are changing along
                    if ( (*arr1)[i]+d1 > p1[axis1] && (*arr1)[i]-d1 < p2[axis1] ) {
                        ret_vec.push_back(i);
                    }
                }
            }

        }   //end for

        return ret_vec;

    }

    [[nodiscard]] inline auto get_middle_inds() const noexcept {
        const double middle_z = ( maxs.z() + mins.z() ) /2;

        std::vector<unsigned long> ret_vec;
        ret_vec.reserve( no_points_unif.x() * no_points_unif.y() );
        for (unsigned long i = 0; i < z.size(); i++) {
            if ( z[i] <= middle_z && z[i]+dz > middle_z ) {
                ret_vec.push_back(i);
            }
        }
        return ret_vec;
    }

    [[nodiscard]] inline auto get_some_x_inds(const double xp) const noexcept {
        std::vector<unsigned long> ret_vec;
        ret_vec.reserve( no_points_unif.x() * no_points_unif.y() );
        for (unsigned long i = 0; i < x.size(); i++) {
            if ( x[i] <= xp && x[i]+dx > xp ) {
                ret_vec.push_back(i);
            }
        }
        return ret_vec;
    }

    void create_no_points_unif() {
        no_points_unif = round( (maxs - mins) / vec3(dx, dy, dz) );
    }

    [[nodiscard]] auto convert_indices_unif(const vec3 &inds) const {
        return static_cast<unsigned long>( inds.x() + inds.y()*no_points_unif.x() + inds.z()*no_points_unif.x()*no_points_unif.y() );
    }

    [[nodiscard]] inline vec3 get_ind_unif(vec3 p) const noexcept {
        p -= mins;
        return {floor(p.x()/dx), floor(p.y()/dy),floor(p.z()/dz)};
    }


    [[nodiscard]] inline bool has_left(const unsigned ind) const noexcept {return r[ind].left != -1;} //returns true if it has a left
    [[nodiscard]] inline bool has_right(const unsigned ind) const noexcept {return r[ind].right != -1;}
    [[nodiscard]] inline bool has_down(const unsigned ind) const noexcept {return r[ind].down != -1;}
    [[nodiscard]] inline bool has_up(const unsigned ind) const noexcept {return r[ind].up != -1;}
    [[nodiscard]] inline bool has_front(const unsigned ind) const noexcept {return r[ind].front != -1;}
    [[nodiscard]] inline bool has_back(const unsigned ind) const noexcept {return r[ind].back != -1;}

    [[nodiscard]] inline bool has_2left(const unsigned ind) const noexcept {
        if (r[ind].left == -1) return false;
        return r[r[ind].left].left != -1;
    }
    [[nodiscard]] inline bool has_2right(const unsigned ind) const noexcept {
        if (r[ind].right == -1) return false;
        return r[r[ind].right].right != -1;
    }
    [[nodiscard]] inline bool has_2down(const unsigned ind) const noexcept {
        if (r[ind].down == -1) return false;
        return r[r[ind].down].down != -1;
    }
    [[nodiscard]] inline bool has_2up(const unsigned ind) const noexcept {
        if (r[ind].up == -1) return false;
        return r[r[ind].up].up != -1;
    }
    [[nodiscard]] inline bool has_2front(const unsigned ind) const noexcept {
        if (r[ind].front == -1) return false;
        return r[r[ind].front].front != -1;
    }
    [[nodiscard]] inline bool has_2back(const unsigned ind) const noexcept {
        if (r[ind].back == -1) return false;
        return r[r[ind].back].back != -1;
    }

    [[nodiscard]] inline bool is_boundary(const unsigned ind) const noexcept {
        //do smarter
        // - check 3by3by region first
        // - then check above and below these points
        //or maybe just list out all the indices
        // - make should to check coordinate axis first before small diagonals before big digaonals
        //only in debug do we check that ind itself is not a boundary (i.e. ind!=-1)
#ifndef NDEBUG
        if (ind==-1) {
            std::cerr << "trying to find if boundary of ind=-1. This should never happen\n";
            return true;
        }
#endif
        if (r[ind].left== -1) {return true;}
        if (r[ind].right== -1) {return true;}
        if (r[ind].down== -1) {return true;}
        if (r[ind].up== -1) {return true;}
        if (r[ind].front== -1) {return true;}
        if (r[ind].back== -1) {return true;}

        if ( r[r[ind].left].up== -1) {return true;}
        if ( r[r[ind].left].down== -1) {return true;}
        if ( r[r[ind].left].front== -1) {return true;}
        if ( r[r[ind].left].back== -1) {return true;}

        if ( r[r[ind].right].up== -1) {return true;}
        if ( r[r[ind].right].down== -1) {return true;}
        if ( r[r[ind].right].front== -1) {return true;}
        if ( r[r[ind].right].back== -1) {return true;}

        if ( r[r[r[ind].left].up].front == -1 ) {return true;}
        if ( r[r[r[ind].left].up].back == -1 ) {return true;}
        if ( r[r[r[ind].left].down].front == -1 ) {return true;}
        if ( r[r[r[ind].left].down].back == -1 ) {return true;}

        if ( r[r[r[ind].right].up].front == -1 ) {return true;}
        if ( r[r[r[ind].right].up].back == -1 ) {return true;}
        if ( r[r[r[ind].right].down].front == -1 ) {return true;}
        if ( r[r[r[ind].right].down].back == -1 ) {return true;}

        if ( r[r[ind].up].front== -1) {return true;}
        if ( r[r[ind].up].back== -1) {return true;}
        if ( r[r[ind].down].front== -1) {return true;}
        if ( r[r[ind].down].back== -1) {return true;}

        /*const bool ends = !has_left(ind) || !has_right(ind) || !has_down(ind) || !has_up(ind) || !has_front(ind) || !has_back(ind);
        if (ends) {
            return true;
        }

        const auto ind_d1 = r[ind].left;
        const auto ind_d2 = r[ind].right;

        const std::array<int, 4> inds_ddiag = { r[r[ind].left].up, r[r[ind].left].down, r[r[ind].right].up, r[r[ind].right].down };
        for (const auto i : inds_ddiag) {
            const bool ddiags = !has_front(i) || !has_back(i);
            if (ddiags) {
                return true;
            }
        }*/

        //return end_x || end_y || end_z;
        return false;
    }

    [[nodiscard]] inline bool is_inside_boundary(const unsigned ind) const noexcept {
        const auto end_x = !has_left(ind) && !has_right(ind);
        const auto end_y = !has_down(ind) && !has_up(ind);
        const auto end_z = !has_front(ind) && !has_back(ind);
        return end_x && end_y && end_z;
    }


    void DEBUG_write_boundary_points() const {
        std::ofstream output("../DEBUG/boundary_points.txt");
        if (output.is_open()) {
            const auto inds = get_middle_inds();
            for (const auto ind : inds) {
                output << x[ind] << " " << y[ind] << " " << is_boundary(ind) << "\n";
            }

        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }

    void DEBUG_write_boundary_points_at_x(const double xp) const {
        std::ofstream output("../DEBUG/boundary_points.txt");
        if (output.is_open()) {
            const auto inds = get_some_x_inds(xp);
            for (const auto ind : inds) {
                output << y[ind] << " " << z[ind] << " " << is_boundary(ind) << "\n";
            }

        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }



};

#endif //CODE_GRID_HPP
