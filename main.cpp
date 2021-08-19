#define EIGEN_NO_AUTOMATIC_RESIZING

#include <iostream>
#include "Rigid_body/body.hpp"
#include <string>

#include "MyMath/ops.hpp"
#include "Examples/Rigid_body_ex.hpp"
#include "Examples/derivatives.hpp"
#include "MyMath/finite_difference.hpp"

#include "MyMath/big_vec.hpp"
#include "MyMath/calc.hpp"

//#include "Examples/big_vec_derivs.h"
//#include "Examples/calc_derivs.hpp"

#include "MyMath/big_matrix.hpp"
#include "Fluid_flow/make_vecs.hpp"
#include "Fluid_flow/make_mats.hpp"

#include "Fluid_flow/create_flow.hpp"

int main() {
    //big_veg_derivs_ex();
    //calc_derivs_ex();

    //big_matrix<2,2,2> mat;
    //mat.add_elm(0,0,0,1,1,1, 2);
    //std::cout << mat << "\n";

    solve_flow();


    /*constexpr double dx = 0.1;
    constexpr double dy = 0.1;
    constexpr double dz = 0.1;
    constexpr double dt = 0.1;
    constexpr double Re = 10;

    big_vec<10,10,10,vec3> b(dx, dy, dz);
    big_vec<10,10,10,vec3> v_n(dx, dy, dz);
    big_vec<10,10,10,vec3> v_n1(dx, dy, dz);
    big_vec<10,10,10,double> p(dx, dy, dz);
    big_vec<10,10,10,double> s(dx, dy, dz);

    big_matrix<10,10,10> Q;
    big_matrix<10,10,10> A;

    std::cout << "Running program useless atm\n";
    make_b_first(b, 100, 0.1, v_n, p);
    make_b(b, 100, 0.1, v_n, v_n1, p);

    make_s_first(s, 100, 0.1, v_n);
    make_s(s, 100, 0.1, v_n, v_n1);

    make_Q(Q, p);
    make_A(A, v_n,dt, Re);*/

    return 0;
}
