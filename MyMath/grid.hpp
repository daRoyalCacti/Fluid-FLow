//
// Created by jacob on 11/9/21.
//

#ifndef CODE_GRID_HPP
#define CODE_GRID_HPP

#include <vector>
#include <limits>
#include <set>
#include <fstream>
#include "../MyMath/vec3.hpp"
#include "../MyMath/dist_to_plane.hpp"

//stores the neighbours of a grid point
// -1 if no neighbour
struct grid_relation {
    int left{}, right{}, down{}, up{}, front{}, back{};

    [[nodiscard]] bool is_edge() const {
        const unsigned x = left==-1 || right==-1;
        const unsigned y = down==-1 || up==-1;
        const unsigned z = front==-1 || back==-1;

        return (x+y+z) > 1;
    }
};

struct axes {
    vec3 x,y,z;
    template <unsigned i>
    [[nodiscard]] vec3& get() {
        if constexpr (i==0) {
            return x;
        } else if constexpr (i==1) {
            return y;
        } else if constexpr( i==2) {
            return z;
        } else {
            std::cerr << "trying to return an axis that doesn't exist\n";
            return z;
        }
    }

    template <unsigned i>
    [[nodiscard]] vec3 get() const {
        if constexpr (i==0) {
            return x;
        } else if constexpr (i==1) {
            return y;
        } else if constexpr( i==2) {
            return z;
        } else {
            std::cerr << "trying to return an axis that doesn't exist\n";
            return z;
        }
    }
};

//grids must be axis aligned
// grid[i] corresponds to the interval x[i]+dx, y[i]+dy, z[i]+dz
struct grid {
    std::vector<double> x, y, z;
    std::vector<double> plot_x, plot_y, plot_z;
    double dx, dy, dz;
    std::vector<grid_relation> r{};
    //vec3 mins, maxs;    //min = (minx, miny, minz)
    vec3 edge1, edge2, /*edge3,*/ edge4, edge5, /*edge6,*/ edge7, edge8;
    vec3 middle;
    vec3 no_points_unif;    //the total number of points across each axis assuming the grid is boring

    axes axis{ vec3(1,0,0), vec3(0,1,0), vec3(0,0,1) };


    vec3 operator[](unsigned i) const noexcept { return { x[i], y[i], z[i] }; }

    grid() : dx{}, dy{}, dz{} {};
    grid(const double dx_, const double dy_, const double dz_, const double minx, const double miny, const double minz, const double maxx, const double maxy, const double maxz)
        : dx(dx_), dy(dy_), dz(dz_), //mins(minx, miny, minz), maxs(maxx, maxy, maxz) {
        middle( (minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2 ),
        edge1(minx, miny, minz), edge2(maxx, miny, minz), /*edge3(minx, maxy, minz),*/ edge4(maxx, maxy, minz),
        edge5(minz, miny, maxz), /*edge6(maxx, miny, maxz),*/ edge7(minx, maxy, maxz), edge8(maxx, maxy, maxz),
        no_points_unif(round( ( vec3(maxx-minx, maxy-miny, maxz-minz) ) / vec3(dx, dy, dz) ) ) {}

    grid(const double dx_, const double dy_, const double dz_, const vec3& minv, const vec3& maxv)
        : grid(dx_, dy_, dz_, minv.x(), minv.y(), minv.z(), maxv.x(), maxv.y(), maxv.z()) {}

    grid(const std::vector<double> &x_, const std::vector<double> &y_, const std::vector<double> &z_, const double dx_, const double dy_, const double dz_,
         const double minx, const double miny, const double minz, const double maxx, const double maxy, const double maxz)
        : x(x_), y(y_), z(z_), dx(dx_), dy(dy_), dz(dz_), plot_x(x_), plot_y(y_), plot_z(z_), //mins(minx, miny, minz), maxs(maxx, maxy, maxz) {
        middle( (minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2 ),
        edge1(minx, miny, minz), edge2(maxx, miny, minz), /*edge3(minx, maxy, minz),*/ edge4(maxx, maxy, minz),
        edge5(minx, miny, maxz), /*edge6(maxx, miny, maxz),*/ edge7(minx, maxy, maxz), edge8(maxx, maxy, maxz),
        no_points_unif(round( ( vec3(maxx-minx, maxy-miny, maxz-minz) ) / vec3(dx, dy, dz) ) ) {
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

    void set_plotting_points() {
        plot_x = x;
        plot_y = y;
        plot_z = z;
    }

    [[nodiscard]] vec3 get_plot_pos(const unsigned ind) const {
        return {plot_x[ind], plot_y[ind], plot_z[ind]};
    }

    void update_pos(const vec3& vel, const vec3& c_o_m, const vec3& omega, const double dt) noexcept {
        // const vec3 rot_angle_vec = model.w*dt;
        //rotate(pos_cm_old, rot_angle_vec, x, rot_angle_vec.length()) + model.v*dt;}



        const auto rot_angle_vec = omega*dt;
        const auto trans_vec = vel*dt;

        for (unsigned i = 0; i < x.size(); i++) {
            const auto rot = rotate(c_o_m, rot_angle_vec, vec3(x[i], y[i], z[i]), rot_angle_vec.length() );
            x[i] = rot.x()+trans_vec.x();
            y[i] = rot.y()+trans_vec.y();
            z[i] = rot.z()+trans_vec.z();
            plot_x[i] += trans_vec.x();
            plot_y[i] += trans_vec.y();
            plot_z[i] += trans_vec.z();
        }

        edge1 = rotate(c_o_m, rot_angle_vec, edge1, rot_angle_vec.length() ) + trans_vec;
        edge2 = rotate(c_o_m, rot_angle_vec, edge2, rot_angle_vec.length() ) + trans_vec;
        //edge3 = rotate(c_o_m, rot_angle_vec, edge3, rot_angle_vec.length() );
        edge4 = rotate(c_o_m, rot_angle_vec, edge4, rot_angle_vec.length() ) + trans_vec;
        edge5 = rotate(c_o_m, rot_angle_vec, edge5, rot_angle_vec.length() ) + trans_vec;
        //edge6 = rotate(c_o_m, rot_angle_vec, edge6, rot_angle_vec.length() );
        edge7 = rotate(c_o_m, rot_angle_vec, edge7, rot_angle_vec.length() ) + trans_vec;
        edge8 = rotate(c_o_m, rot_angle_vec, edge8, rot_angle_vec.length() ) + trans_vec;
        middle = rotate(c_o_m, rot_angle_vec, middle, rot_angle_vec.length() ) + trans_vec;

        axis.x = rotate(vec3(0), -rot_angle_vec, axis.x, rot_angle_vec.length() );
        axis.y = rotate(vec3(0), -rot_angle_vec, axis.y, rot_angle_vec.length() );
        axis.z = rotate(vec3(0), -rot_angle_vec, axis.z, rot_angle_vec.length() );



        /*for (auto & x_ : x) {
            x_ += x_off;
        }
        for (auto & y_ : y) {
            y_ += y_off;
        }
        for (auto & z_ : z) {
            z_ += z_off;
        }
        const vec3 off(x_off, y_off, z_off);
        mins+=off;  //might need to store positions of corners - should only actualy need to store 4 of them
        maxs+=off;  //needed for off_wall*/
    }

    //check to see if an index is away from the walls
    [[nodiscard]] bool off_walls(const unsigned i) const {
        /*const bool away_x = (x[i] > mins.x()+2*dx) && (x[i]<maxs.x()-2*dx);
        const bool away_y = (y[i] > mins.y()+2*dy) && (y[i]<maxs.y()-2*dy);
        const bool away_z = (z[i] > mins.z()+2*dz) && (z[i]<maxs.z()-2*dz);*/
        const vec3 v = {x[i], y[i], z[i]};
        //the distances to all the walls
        const auto dist_left = dist_to_plane( v, edge1, edge5, edge7 );
        const auto dist_right = dist_to_plane( v, edge2, edge4, edge8 );
        const auto dist_down = dist_to_plane( v, edge1, edge2, edge5 );
        const auto dist_up = dist_to_plane( v, edge4, edge8, edge7 );
        const auto dist_front = dist_to_plane( v, edge1, edge2, edge4 );
        const auto dist_back = dist_to_plane( v, edge5, edge7, edge8 );

        //std::cerr << i << "\t" << v << "\n";
        //std::cerr << "\t" << dist_left << " " << dist_right << " " << dist_down << " " << dist_up << " " << dist_front << " " << dist_back << "\n";

        const bool away_x = (dist_left > 2*dx) && (dist_right > 2*dx);
        const bool away_y = (dist_down > 2*dy) && (dist_up > 2*dx);
        const bool away_z = (dist_front > 2*dz) && (dist_back > 2*dz);


        return away_x && away_y && away_z;
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

    [[maybe_unused]] [[nodiscard]] inline auto get_inds(vec3 p1, vec3 p2) const noexcept {
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
        //const double middle_z = ( maxs.z() + mins.z() ) /2;
        const auto middle_z = middle.z();

        std::vector<unsigned long> ret_vec;
        ret_vec.reserve( static_cast<unsigned long>(no_points_unif.x() * no_points_unif.y()) );
        for (unsigned long i = 0; i < z.size(); i++) {
            //NOPE JUST GET INDS INITIALLY AND USE THAT FOR ALL WRITING
            if ( z[i] <= middle_z && z[i] >= middle_z-dz ) {    //TODO : update this to used the max of dx,dy,dz (when the mesh rotates, dx could be in the z direction)
                ret_vec.push_back(i);
            }
        }
        return ret_vec;
    }

    [[nodiscard]] inline auto get_some_x_inds(const double xp) const noexcept {
        std::vector<unsigned long> ret_vec;
        ret_vec.reserve( static_cast<unsigned long>(no_points_unif.x() * no_points_unif.y()) );
        for (unsigned long i = 0; i < x.size(); i++) {
            if ( x[i] <= xp && x[i]+dx > xp ) {
                ret_vec.push_back(i);
            }
        }
        return ret_vec;
    }

    /*inline void create_no_points_unif() {
        no_points_unif = round( (maxs - mins) / vec3(dx, dy, dz) );
    }*/

    [[nodiscard]] auto convert_indices_unif(const vec3 &inds) const {
        return static_cast<unsigned long>( inds.x() + inds.y()*no_points_unif.x() + inds.z()*no_points_unif.x()*no_points_unif.y() );
    }

    /*[[nodiscard]] inline vec3 get_ind_unif(vec3 p) const noexcept {
        p -= mins;
        return {floor(p.x()/dx), floor(p.y()/dy),floor(p.z()/dz)};
    }*/


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


        /*if ( r[r[ind].left].up== -1) {return true;}
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
        if ( r[r[ind].down].back== -1) {return true;}*/


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


    void DEBUG_write_boundary_points(const bool print_off_wall = false) const {
        std::ofstream output("../DEBUG/boundary_points.txt");
        if (output.is_open()) {
            const auto inds = get_middle_inds();
            for (const auto ind : inds) {
                if (!print_off_wall) {
                    output << x[ind] << " " << y[ind] << " " << is_boundary(ind) << "\n";
                } else {
                    output << x[ind] << " " << y[ind] << " " << (is_boundary(ind) && off_walls(ind)) << "\n";
                }
            }

        } else {
            std::cerr << "failed to open file\n";
        }

        output.close();
    }

    [[maybe_unused]] void DEBUG_write_boundary_points_at_x(const double xp) const {
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



    [[nodiscard]] bool can_move(const unsigned ind, const int x_m, const int y_m, const int z_m) const noexcept {
        auto curr_ind = ind;
#ifndef NDEBUG
        if (curr_ind == -1) {
            std::cerr << "index is -1. This should never happen.\n";
        }
        if (curr_ind >= size()) {
            std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
        }
#endif

        if (x_m > 0) {
            for (unsigned i = 0; i < x_m; i++) {
                curr_ind = r[curr_ind].right;
                if (curr_ind == -1) {
                    return false;
                }
            }
        } else {
            for (int i = x_m; i < 0; i++) {
                curr_ind = r[curr_ind].left;
                if (curr_ind == -1) {
                    return false;
                }
            }
        }
        if (y_m > 0) {
            for (unsigned j = 0; j < y_m; j++) {
                curr_ind = r[curr_ind].up;
                if (curr_ind == -1) {
                    return false;
                }
            }
        } else {
            for (int j = y_m; j < 0; j++) {
                curr_ind = r[curr_ind].down;
                if (curr_ind == -1) {
                    return false;
                }
            }
        }
        if (z_m > 0) {
            for (unsigned k = 0; k < z_m; k++) {
                curr_ind = r[curr_ind].back;
                if (curr_ind == -1) {
                    return false;
                }
            }
        } else {
            for (int k = z_m; k < 0; k++) {
                curr_ind = r[curr_ind].front;
                if (curr_ind == -1) {
                    return false;
                }
            }
        }

        return true;
    }

    [[nodiscard]] auto can_move(const unsigned ind, const vec3&  v) const noexcept {
        return can_move(ind, static_cast<int>(v.x()), static_cast<int>(v.y()), static_cast<int>(v.z()) );
    }


    [[nodiscard]] unsigned get_move_ind(const unsigned ind, const int x_m, const int y_m, const int z_m) const noexcept {
        auto curr_ind = ind;
#ifndef NDEBUG
        if (curr_ind == -1) {
            std::cerr << "index is -1. This should never happen.\n";
        }
        if (curr_ind >= size()) {
            std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
        }

        if (!can_move(ind, x_m, y_m, z_m)) {
            std::cerr << "trying to move to a point that cannot be moved to\n";
        }
#endif

        if (x_m > 0) {
            for (unsigned i = 0; i < x_m; i++) {
                curr_ind = r[curr_ind].right;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\tx_m move failed. current x_m : " << i+1 << "/" << x_m << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tx_m move failed. current x_m : " << i+1 << "/" << x_m << "\n";
                }
#endif
            }
        } else {
            for (int i = x_m; i < 0; i++) {
                curr_ind = r[curr_ind].left;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\tx_m move failed. current x_m : " << i+1 << "/" << x_m << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tx_m move failed. current x_m : " << i+1 << "/" << x_m << "\n";
                }
#endif
            }
        }
        if (y_m > 0) {
            for (unsigned j = 0; j < y_m; j++) {
                curr_ind = r[curr_ind].up;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\ty_m move failed. current y_m : " << j+1 << "/" << y_m << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\ty_m move failed. current y_m : " << j+1 << "/" << y_m << "\n";
                }
#endif
            }
        } else {
            for (int j = y_m; j < 0; j++) {
                curr_ind = r[curr_ind].down;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\ty_m move failed. current y_m : " << j+1 << "/" << y_m << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\ty_m move failed. current y_m : " << j+1 << "/" << y_m << "\n";
                }
#endif
            }
        }
        if (z_m > 0) {
            for (unsigned k = 0; k < z_m; k++) {
                curr_ind = r[curr_ind].back;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\tz_m move failed. current z_m : " << k+1 << "/" << z_m << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tz_m move failed. current z_m : " << k+1 << "/" << z_m << "\n";
                }
#endif
            }
        } else {
            for (int k = z_m; k < 0; k++) {
                curr_ind = r[curr_ind].front;
#ifndef NDEBUG
                if (curr_ind == -1) {
                    std::cerr << "current index is -1. This should never happen.\n";
                    std::cerr << "\tz_m move failed. current z_m : " << k+1 << "/" << z_m << "\n";
                }
                if (curr_ind >= size()) {
                    std::cerr << "index is outside the range of the grid\n\tcurr_ind=" << curr_ind << " size=" << size() << "\n";
                    std::cerr << "\tz_m move failed. current z_m : " << k+1 << "/" << z_m << "\n";
                }
#endif
            }
        }

        return curr_ind;
    }

    [[nodiscard]] auto get_move_ind(const unsigned ind, const vec3&  v) const noexcept {
        return get_move_ind(ind, static_cast<int>(v.x()), static_cast<int>(v.y()), static_cast<int>(v.z()) );
    }


};

#endif //CODE_GRID_HPP
