#define EIGEN_NO_AUTOMATIC_RESIZING

/*
#include <iostream>
#include "Rigid_body/body.hpp"
#include <string>

#include "MyMath/ops.hpp"
#include "Examples/Rigid_body_ex.hpp"
#include "Examples/derivatives.hpp"
#include "MyMath/finite_difference.hpp"

#include "MyMath/big_vec.hpp"
#include "MyMath/calc.hpp"
*/

#include "Examples/big_vec_derivs.h"
#include "Examples/calc_derivs.hpp"
/*
#include "MyMath/big_matrix.hpp"
#include "Fluid_flow/make_vecs.hpp"
#include "Fluid_flow/make_mats.hpp"

#include "Fluid_flow/create_flow.hpp"
*/

int main() {
    big_veg_derivs_ex();
    calc_derivs_ex();

    //solve_flow();



    return 0;
}
