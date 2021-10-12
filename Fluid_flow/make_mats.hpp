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
void make_Q(big_matrix &Q, const big_vec_d &p, const boundary_conditions &bc) noexcept {
    const auto dx = p.g->dx;
    const auto dy = p.g->dy;
    const auto dz = p.g->dz;

    const auto &norms = bc.norms;


//#pragma omp parallel for
    //shared(Q, p, dxdx, dydy, dzdz) default(none)
    for (unsigned i = 0; i < p.size(); i++) {

        if (p.is_boundary(i)) {
            if (p.g->off_walls(i)) {
                const auto &n =  bc.norms.normal(i);
                double mid = 0;

                //x
                if (p.can_move(i, -1,0,0) && p.can_move(i, 1,0,0)) {    //central difference
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), n.x()*1/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, -1,0,0), -n.x()*1/(2*dx) );
                } else if (p.can_move(i, 2,0,0)) {  //forward
                    mid += -n.x()*3/(2*dx);
                    Q.add_elm(i, p.get_move_ind(i, 2,0,0), -n.x()*1/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), n.x()*4/(2*dx) );
                } else {    //backward
                    mid += n.x()*3/(2*dx);
                    Q.add_elm(i, p.get_move_ind(i, -1,0,0), -n.x()*4/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, -2,0,0), n.x()*1/(2*dx) );
                }

                //y
                if (p.can_move(i, 0,-1,0) && p.can_move(i, 0,1,0)) {    //central difference
                    Q.add_elm(i, p.get_move_ind(i, 0,1,0), n.y()*1/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,-1,0), -n.y()*1/(2*dy) );
                } else if (p.can_move(i, 0,2,0)) {  //forward
                    mid += -n.y()*3/(2*dy);
                    Q.add_elm(i, p.get_move_ind(i, 0,2,0), -n.y()*1/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,1,0), n.y()*4/(2*dy) );
                } else {    //backward
                    mid += n.y()*3/(2*dy);
                    Q.add_elm(i, p.get_move_ind(i, 0,-1,0), -n.y()*4/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,-2,0), n.y()*1/(2*dy) );
                }

                //z
                if (p.can_move(i, 0,0,-1) && p.can_move(i, 0,0,1)) {    //central difference
                    Q.add_elm(i, p.get_move_ind(i, 0,0,1), n.z()*1/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-1), -n.z()*1/(2*dz) );
                } else if (p.can_move(i, 0,0,2)) {  //forward
                    mid += -n.z()*3/(2*dz);
                    Q.add_elm(i, p.get_move_ind(i, 0,0,2), -n.z()*1/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,1), n.z()*4/(2*dz) );
                } else {    //backward
                    mid += n.z()*3/(2*dz);
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-1), -n.z()*4/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-2), n.z()*1/(2*dz) );
                }

                /*if (p.can_move(i, -1,0,0) && p.can_move(i, 1,0,0)) {    //central difference
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), n.x()*1/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, -1,0,0), -n.x()*1/(2*dx) );
                } else if (p.can_move(i, 1,0,0)) {  //forward
                    mid += -n.x()/dx;
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), n.x()/dx );
                } else {    //backward
                    mid += n.x()/dx;
                    Q.add_elm(i, p.get_move_ind(i, -1,0,0), -n.x()/dx );
                }

                //y
                if (p.can_move(i, 0,-1,0) && p.can_move(i, 0,1,0)) {    //central difference
                    Q.add_elm(i, p.get_move_ind(i, 0,1,0), n.y()*1/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,-1,0), -n.y()*1/(2*dy) );
                } else if (p.can_move(i, 0,1,0)) {  //forward
                    mid += -n.y()/dy;
                    Q.add_elm(i, p.get_move_ind(i, 0,1,0), n.y()/dy );
                } else {    //backward
                    mid += n.y()*3/(2*dz);
                    Q.add_elm(i, p.get_move_ind(i, 0,-1,0), -n.y()/dy );
                }

                //z
                if (p.can_move(i, 0,0,-1) && p.can_move(i, 0,0,1)) {    //central difference
                    Q.add_elm(i, p.get_move_ind(i, 0,0,1), n.z()*1/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-1), -n.z()*1/(2*dz) );
                } else if (p.can_move(i, 0,0,1)) {  //forward
                    mid += -n.z()/dz;
                    Q.add_elm(i, p.get_move_ind(i, 0,0,1), n.z()/dz );
                } else {    //backward
                    mid += n.z()/dz;
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-1), -n.z()/dz );
                }*/


                Q.add_elm(i,i, mid);
            } else {    //on walls
                Q.add_elm(i, i,1);
                //d^2p/dn^2=0
                /*const auto &n =  bc.norms.normal(i);
#ifndef NDEBUG
                if (n != vec3(1,0,0) && n != vec3(-1,0,0) &&
                n!= vec3(0,1,0) && n!= vec3(0,-1,0) &&
                n!= vec3(0,0,1) && n != vec3(0,0,-1)) {
                    std::cerr << "wall normal set incorrectly\n";
                }
#endif
//first order derivative
//A.add_elm(i,i, 1);
//A.add_elm(i,v.get_move_ind(i,n), -1);



                unsigned dir;
#ifndef NDEBUG
                bool set = false;
#endif
                for (unsigned j = 0; j < 3; j++) {
                    if (n[j] != 0 ) {
                        dir = j;
#ifndef NDEBUG
                        set = true;
#endif
                    }
                }
#ifndef NDEBUG
                if (!set) {
                    std::cerr << "wall normal is the zero vector\n";
                }
#endif
                double hh;
                if (dir == 0) {
                    hh = dx*dx;
                } else if (dir == 1) {
                    hh = dy*dy;
                } else {
                    hh = dz*dz;
                }

                Q.add_elm(i,i, 2/hh);
                Q.add_elm(i, p.get_move_ind(i, n), -5/hh);
                Q.add_elm(i, p.get_move_ind(i, 2*n), 4/hh);
                Q.add_elm(i, p.get_move_ind(i, 3*n), -1/hh);
                */
            }
        } else {

            const auto dxdx = p.dx(i)*p.dx(i);
            const auto dydy = p.dy(i)*p.dy(i);
            const auto dzdz = p.dz(i)*p.dz(i);


            //x-axis
            Q.add_elm(i, p.get_move_ind(i,1,0,0), 1 / dxdx);
            Q.add_elm(i, p.get_move_ind(i,-1,0,0), 1 / dxdx);


            //y-axis
            Q.add_elm(i, p.get_move_ind(i,0,1,0), 1 / dydy);
            Q.add_elm(i, p.get_move_ind(i,0,-1,0), 1 / dydy);


            //z-axis
            Q.add_elm(i, p.get_move_ind(i,0,0,1), 1 / dzdz);
            Q.add_elm(i, p.get_move_ind(i,0,0,-1), 1 / dzdz);


            Q.add_elm(i, i, -2/dxdx - 2/dydy - 2/dzdz);

        }

    }
}



void make_A(big_matrix &A,  const big_vec_v &v, const double dt, const double Re, const boundary_conditions &bc) noexcept {
    const auto dx = v.g->dx;
    const auto dy = v.g->dy;
    const auto dz = v.g->dz;

    //#pragma omp parallel for
    //shared(A, v, dt, Re, Rdxdx, Rdydy, Rdzdz) default(none)
    for (unsigned i = 0; i < v.g->size(); i++) {
        if (v.is_boundary(i)) {
            if (v.g->off_walls(i)) {
                A.add_elm(i,i, 1);
            } else {
                const auto &n =  bc.norms.normal(i);
                double mid = 0;
                //x
                if (v.can_move(i, -1,0,0) && v.can_move(i, 1,0,0)) {    //central difference
                    A.add_elm(i, v.get_move_ind(i, 1,0,0), n.x()*1/(2*dx) );
                    A.add_elm(i, v.get_move_ind(i, -1,0,0), -n.x()*1/(2*dx) );
                } else if (v.can_move(i, 2,0,0)) {  //forward
                    mid += -n.x()*3/(2*dx);
                    A.add_elm(i, v.get_move_ind(i, 2,0,0), -n.x()*1/(2*dx) );
                    A.add_elm(i, v.get_move_ind(i, 1,0,0), n.x()*4/(2*dx) );
                } else {    //backward
                    mid += n.x()*3/(2*dx);
                    A.add_elm(i, v.get_move_ind(i, -1,0,0), -n.x()*4/(2*dx) );
                    A.add_elm(i, v.get_move_ind(i, -2,0,0), n.x()*1/(2*dx) );
                }

                //y
                if (v.can_move(i, 0,-1,0) && v.can_move(i, 0,1,0)) {    //central difference
                    A.add_elm(i, v.get_move_ind(i, 0,1,0), n.y()*1/(2*dy) );
                    A.add_elm(i, v.get_move_ind(i, 0,-1,0), -n.y()*1/(2*dy) );
                } else if (v.can_move(i, 0,2,0)) {  //forward
                    mid += -n.y()*3/(2*dy);
                    A.add_elm(i, v.get_move_ind(i, 0,2,0), -n.y()*1/(2*dy) );
                    A.add_elm(i, v.get_move_ind(i, 0,1,0), n.y()*4/(2*dy) );
                } else {    //backward
                    mid += n.y()*3/(2*dy);
                    A.add_elm(i, v.get_move_ind(i, 0,-1,0), -n.y()*4/(2*dy) );
                    A.add_elm(i, v.get_move_ind(i, 0,-2,0), n.y()*1/(2*dy) );
                }

                //z
                if (v.can_move(i, 0,0,-1) && v.can_move(i, 0,0,1)) {    //central difference
                    A.add_elm(i, v.get_move_ind(i, 0,0,1), n.z()*1/(2*dz) );
                    A.add_elm(i, v.get_move_ind(i, 0,0,-1), -n.z()*1/(2*dz) );
                } else if (v.can_move(i, 0,0,2)) {  //forward
                    mid += -n.z()*3/(2*dz);
                    A.add_elm(i, v.get_move_ind(i, 0,0,2), -n.z()*1/(2*dz) );
                    A.add_elm(i, v.get_move_ind(i, 0,0,1), n.z()*4/(2*dz) );
                } else {    //backward
                    mid += n.z()*3/(2*dz);
                    A.add_elm(i, v.get_move_ind(i, 0,0,-1), -n.z()*4/(2*dz) );
                    A.add_elm(i, v.get_move_ind(i, 0,0,-2), n.z()*1/(2*dz) );
                }

                /*if (v.can_move(i, -1,0,0) && v.can_move(i, 1,0,0)) {    //central difference
                    A.add_elm(i, v.get_move_ind(i, 1,0,0), n.x()*1/(2*dx) );
                    A.add_elm(i, v.get_move_ind(i, -1,0,0), -n.x()*1/(2*dx) );
                } else if (v.can_move(i, 1,0,0)) {  //forward
                    mid += -n.x()/dx;
                    A.add_elm(i, v.get_move_ind(i, 1,0,0), n.x()/dx );
                } else {    //backward
                    mid += n.x()/dx;
                    A.add_elm(i, v.get_move_ind(i, -1,0,0), -n.x()/dx );
                }

                //y
                if (v.can_move(i, 0,-1,0) && v.can_move(i, 0,1,0)) {    //central difference
                    A.add_elm(i, v.get_move_ind(i, 0,1,0), n.y()*1/(2*dy) );
                    A.add_elm(i, v.get_move_ind(i, 0,-1,0), -n.y()*1/(2*dy) );
                } else if (v.can_move(i, 0,1,0)) {  //forward
                    mid += -n.y()/dy;
                    A.add_elm(i, v.get_move_ind(i, 0,1,0), n.y()/dy );
                } else {    //backward
                    mid += n.y()*3/(2*dz);
                    A.add_elm(i, v.get_move_ind(i, 0,-1,0), -n.y()/dy );
                }

                //z
                if (v.can_move(i, 0,0,-1) && v.can_move(i, 0,0,1)) {    //central difference
                    A.add_elm(i, v.get_move_ind(i, 0,0,1), n.z()*1/(2*dz) );
                    A.add_elm(i, v.get_move_ind(i, 0,0,-1), -n.z()*1/(2*dz) );
                } else if (v.can_move(i, 0,0,1)) {  //forward
                    mid += -n.z()/dz;
                    A.add_elm(i, v.get_move_ind(i, 0,0,1), n.z()/dz );
                } else {    //backward
                    mid += n.z()/dz;
                    A.add_elm(i, v.get_move_ind(i, 0,0,-1), -n.z()/dz );
                }*/


                A.add_elm(i,i, mid);

                /*
                const auto &n =  bc.norms.normal(i);
#ifndef NDEBUG
                if (n != vec3(1,0,0) && n != vec3(-1,0,0) &&
                    n!= vec3(0,1,0) && n!= vec3(0,-1,0) &&
                    n!= vec3(0,0,1) && n != vec3(0,0,-1)) {
                    std::cerr << "wall normal set incorrectly\n";
                }
#endif
                //first order derivative
                //A.add_elm(i,i, 1);
                //A.add_elm(i,v.get_move_ind(i,n), -1);



                unsigned dir;
#ifndef NDEBUG
                bool set = false;
#endif
                for (unsigned j = 0; j < 3; j++) {
                    if (n[j] != 0 ) {
                        dir = j;
#ifndef NDEBUG
                        set = true;
#endif
                    }
                }
#ifndef NDEBUG
                if (!set) {
                    std::cerr << "wall normal is the zero vector\n";
                }
#endif
                double h;
                if (dir == 0) {
                    h = v.g->dx;
                } else if (dir == 1) {
                    h = v.g->dy;
                } else {
                    h = v.g-> dz;
                }

                A.add_elm(i,i, 3/(2*h));
                A.add_elm(i, v.get_move_ind(i,n), -2/h);
                A.add_elm(i, v.get_move_ind(i, 2*n), 1/(2*h) );
                 */
            }
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
