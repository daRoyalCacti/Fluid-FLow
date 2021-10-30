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
#ifndef NDEBUG
    if (dx == 0) {
        std::cerr << "Make Q has dx=0\n";
    }
    if (dy == 0) {
        std::cerr << "Make Q has dy=0\n";
    }
    if (dz == 0) {
        std::cerr << "Make Q has dz=0\n";
    }
#endif

    const auto &norms = bc.norms;


    #pragma omp parallel for
    for (unsigned i = 0; i < p.size(); i++) {
        //std::cerr << i << "/" << p.size() << "\n";
       // std::cerr << "\t" << p.is_boundary(i) << "\t" << p.g->off_walls(i) << "\n";

        if (p.is_boundary(i)) {
            if (p.g->off_walls(i)) {
               // std::cerr << "\t\t" << bc.norms.contains(i) << "\n";
                const auto &n =  bc.norms.normal(i);
                double mid = 0;

                //x
                if (p.can_move(i, -1,0,0) && p.can_move(i, 1,0,0)) {    //central difference
                    //std::cerr << "\t\tcentral x" << "\n";
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), n.x()*1/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, -1,0,0), -n.x()*1/(2*dx) );
                } else if (p.can_move(i, 2,0,0)) {  //forward
                    //std::cerr << "\t\tforward x" << "\n";
                    mid += -n.x()*3/(2*dx);
                    Q.add_elm(i, p.get_move_ind(i, 2,0,0), -n.x()*1/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, 1,0,0), n.x()*4/(2*dx) );
                } else if (p.can_move(i, -2,0,0)){    //backward
                    //std::cerr << "\t\tbackward x" << "\n";
                    mid += n.x()*3/(2*dx);
                    Q.add_elm(i, p.get_move_ind(i, -1,0,0), -n.x()*4/(2*dx) );
                    Q.add_elm(i, p.get_move_ind(i, -2,0,0), n.x()*1/(2*dx) );
                } else {
                    std::cerr << "throwing for i = " << i << "\n";
                    std::cerr << "\tthis is at (" << vec3(p.g->x[i], p.g->y[i], p.g->z[i]) << ")\n";
                    throw std::runtime_error("cannot compute a derivative in x");
                }
                //std::cerr << "\t\tend x\n";

                //y
                if (p.can_move(i, 0,-1,0) && p.can_move(i, 0,1,0)) {    //central difference
                    //std::cerr << "\t\tcentral y" << "\n";
                    Q.add_elm(i, p.get_move_ind(i, 0,1,0), n.y()*1/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,-1,0), -n.y()*1/(2*dy) );
                } else if (p.can_move(i, 0,2,0)) {  //forward
                    //std::cerr << "\t\tforward y" << "\n";
                    mid += -n.y()*3/(2*dy);
                    Q.add_elm(i, p.get_move_ind(i, 0,2,0), -n.y()*1/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,1,0), n.y()*4/(2*dy) );
                } else if (p.can_move(i, 0,-2,0)){    //backward
                    //std::cerr << "\t\tbackward y" << "\n";
                    mid += n.y()*3/(2*dy);
                    Q.add_elm(i, p.get_move_ind(i, 0,-1,0), -n.y()*4/(2*dy) );
                    Q.add_elm(i, p.get_move_ind(i, 0,-2,0), n.y()*1/(2*dy) );
                }else {
                    std::cerr << "throwing for i = " << i << "\n";
                    std::cerr << "\tthis is at (" << vec3(p.g->x[i], p.g->y[i], p.g->z[i]) << ")\n";
                    throw std::runtime_error("cannot compute a derivative in y");
                }
                //std::cerr << "\t\tend y\n";

                //z
                if (p.can_move(i, 0,0,-1) && p.can_move(i, 0,0,1)) {    //central difference
                    //std::cerr << "\t\tcentral z" << "\n";
                    Q.add_elm(i, p.get_move_ind(i, 0,0,1), n.z()*1/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-1), -n.z()*1/(2*dz) );
                } else if (p.can_move(i, 0,0,2)) {  //forward
                    //std::cerr << "\t\tforward z" << "\n";
                    mid += -n.z()*3/(2*dz);
                    Q.add_elm(i, p.get_move_ind(i, 0,0,2), -n.z()*1/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,1), n.z()*4/(2*dz) );
                } else if (p.can_move(i, 0,0,-2)){    //backward
                    //std::cerr << "\t\tbackward z" << "\n";
                    mid += n.z()*3/(2*dz);
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-1), -n.z()*4/(2*dz) );
                    Q.add_elm(i, p.get_move_ind(i, 0,0,-2), n.z()*1/(2*dz) );
                }else {
                    std::cerr << "throwing for i = " << i << "\n";
                    std::cerr << "\tthis is at (" << vec3(p.g->x[i], p.g->y[i], p.g->z[i]) << ")\n";
                    throw std::runtime_error("cannot compute a derivative in z");
                }
                //std::cerr << "\t\tend z\n";


                Q.add_elm(i,i, mid);
            } else {    //on walls
                Q.add_elm(i, i,1);

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


    #pragma omp parallel for
    for (unsigned i = 0; i < v.g->size(); i++) {

#ifndef NDEBUG
        if (v.dx(i) == 0) {
            std::cerr << "Make A has dx=0\n";
        }
        if (v.dy(i) == 0) {
            std::cerr << "Make A has dy=0\n";
        }
        if (v.dz(i) == 0) {
            std::cerr << "Make A has dz=0\n";
        }
#endif

        //std::cerr << "\t" << i << "/" << v.g->size() << "\n";
        //std::cerr << "\t\t" << v.is_boundary(i) << "\n";
        if (v.is_boundary(i)) {
            if (v.g->off_walls(i)) {
                A.add_elm(i,i, 1);
            } else {
                //std::cerr << "\t\t" << bc.norms.contains(i) << "\n";

                const auto &n =  bc.norms.normal(i);
#ifndef NDEBUG
                if (n != vec3(1,0,0) && n != vec3(-1,0,0) &&
                    n!= vec3(0,1,0) && n!= vec3(0,-1,0) &&
                    n!= vec3(0,0,1) && n != vec3(0,0,-1)) {
                    std::cerr << "wall normal set incorrectly\n";
                }
#endif

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
