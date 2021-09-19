//
// Created by jacob on 13/8/21.
//

#ifndef CODE_MAKE_MATS_HPP
#define CODE_MAKE_MATS_HPP

#include "../MyMath/big_matrix.hpp"
#include "../MyMath/big_vec.hpp"
#include "../MyMath/boundary.hpp"



//p to store dx,dy,dz and BC

//p contains all boundary info
//need normal infor
// - if points are on boundary, Q is set differenctly
// - which boundary matters depends on the direction of the normal vector
void make_Q(big_matrix &Q, const big_vec_d &p, const boundary_normals &norms) noexcept {



//#pragma omp parallel for
    //shared(Q, p, dxdx, dydy, dzdz) default(none)
    for (unsigned i = 0; i < p.size(); i++) {
        double diag = 0.0;

        if (p.is_boundary(i)) {
#ifndef NDEBUG
            if (!norms.contains(i)) {
                std::cerr << "norms does not contain " << i  << "\n";
            }
#endif
            Q.add_elm(i, i,1);
            /*const auto norm = norms.normal(i);
            //picking the direction
            unsigned big_dir = 0;
            if (std::abs(norm.y()) > std::abs(norm[big_dir]) ) {
                big_dir = 1;
            }
            if (std::abs(norm.z()) > std::abs(norm[big_dir]) ) {
                big_dir = 2;
            }

            if (big_dir == 0) { //x direction biggest
                if (p.has_left(i) && p.has_right(i)) {
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), 1);
                } else {
                    Q.add_elm(i, i,1);
                }
            } else if (big_dir == 1) {  //y direction biggest
                if (p.has_down(i) && p.has_up(i)) {
                    Q.add_elm(i, p.get_move_ind(i, 0, 1, 0),1);
                } else {
                    Q.add_elm(i, i,1);
                }
            } else {    //z direction
#ifndef NDEBUG
                if (big_dir > 2) {
                    std::cerr << "the biggest direction cannot be larger than 2\n";
                }
#endif
                if (p.has_front(i) && p.has_back(i)) {
                    Q.add_elm(i, p.get_move_ind(i,0,0,1),1);
                } else {
                    Q.add_elm(i, i,1);
                }
            }*/

        } else {

            const auto dxdx = p.dx(i)*p.dx(i);
            const auto dydy = p.dy(i)*p.dy(i);
            const auto dzdz = p.dz(i)*p.dz(i);


            //x axis checks
            if (!p.has_left(i)) {
                diag += 2 / dxdx;

                Q.add_elm(i, p.get_move_ind(i,3,0,0), -1 / dxdx);
                Q.add_elm(i, p.get_move_ind(i,2,0,0), 4 / dxdx);
                Q.add_elm(i, p.get_move_ind(i,1,0,0), -5 / dxdx);
            } else if (!p.has_right(i)) {
                diag += 2 / dxdx;

                Q.add_elm(i, p.get_move_ind(i,-1,0,0), -5 / dxdx);
                Q.add_elm(i, p.get_move_ind(i,-2,0,0), 4 / dxdx);
                Q.add_elm(i, p.get_move_ind(i,-3,0,0), -1 / dxdx);
            } else {
                diag += -2 / dxdx;

                Q.add_elm(i, p.get_move_ind(i,1,0,0), 1 / dxdx);
                Q.add_elm(i, p.get_move_ind(i,-1,0,0), 1 / dxdx);
            }

            //y axis checks
            if (!p.has_down(i)) {
                diag += 2 / dydy;

                Q.add_elm(i, p.get_move_ind(i,0,3,0), -1 / dydy);
                Q.add_elm(i, p.get_move_ind(i,0,2,0), 4 / dydy);
                Q.add_elm(i, p.get_move_ind(i,0,1,0), -5 / dydy);
            } else if (!p.has_up(i)) {
                diag += 2 / dydy;

                Q.add_elm(i, p.get_move_ind(i,0,-1,0), -5 / dydy);
                Q.add_elm(i, p.get_move_ind(i,0,-2,0), 4 / dydy);
                Q.add_elm(i, p.get_move_ind(i,0,-3,0), -1 / dydy);
            } else {
                diag += -2 / dydy;

                Q.add_elm(i, p.get_move_ind(i,0,1,0), 1 / dydy);
                Q.add_elm(i, p.get_move_ind(i,0,-1,0), 1 / dydy);
            }

            //z axis checks
            if (!p.has_front(i)) {
                diag += 2 / dzdz;

                Q.add_elm(i, p.get_move_ind(i,0,0,3), -1 / dzdz);
                Q.add_elm(i, p.get_move_ind(i,0,0,2), 4 / dzdz);
                Q.add_elm(i, p.get_move_ind(i,0,0,1), -5 / dzdz);
            } else if (!p.has_back(i)) {
                diag += 2 / dzdz;

                Q.add_elm(i, p.get_move_ind(i,0,0,-1), -5 / dzdz);
                Q.add_elm(i, p.get_move_ind(i,0,0,-2), 4 / dzdz);
                Q.add_elm(i, p.get_move_ind(i,0,0,-3), -1 / dzdz);
            } else {
                diag += -2 / dzdz;

                Q.add_elm(i, p.get_move_ind(i,0,0,1), 1 / dzdz);
                Q.add_elm(i, p.get_move_ind(i,0,0,-1), 1 / dzdz);
            }

            Q.add_elm(i, i, diag);

        }




    }
}



void make_A(big_matrix &A,  const big_vec_v &v, const double dt, const double Re) noexcept {


    //#pragma omp parallel for
    //shared(A, v, dt, Re, Rdxdx, Rdydy, Rdzdz) default(none)
    for (unsigned i = 0; i < v.g->size(); i++) {
        if (v.is_boundary(i)) {
            A.add_elm(i, i,  1);
        } else {
            const auto Rdxdx = Re*v.dx(i)*v.dx(i);
            const auto Rdydy = Re*v.dy(i)*v.dy(i);
            const auto Rdzdz = Re*v.dz(i)*v.dz(i);

            A.add_elm(i,  i,  1/dt + 1/Rdxdx + 1/Rdydy + 1/Rdzdz );
            A.add_elm(i,  v.get_move_ind(i, 1,0,0),  -1/(2*Rdxdx));
            A.add_elm(i,  v.get_move_ind(i, -1,0,0),  -1/(2*Rdxdx));
            A.add_elm(i,  v.get_move_ind(i, 0,1,0),  -1/(2*Rdydy));
            A.add_elm(i,  v.get_move_ind(i, 0,-1,0),  -1/(2*Rdydy));
            A.add_elm(i,  v.get_move_ind(i, 0,0,1),  -1/(2*Rdzdz));
            A.add_elm(i,  v.get_move_ind(i, 0,0,-1),  -1/(2*Rdzdz));
        }


    }
}

#endif //CODE_MAKE_MATS_HPP
