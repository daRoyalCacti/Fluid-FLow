#define EIGEN_NO_AUTOMATIC_RESIZING

//#include "Examples/big_vec_derivs.h"
//#include "Examples/calc_derivs.hpp"


#include "Fluid_flow/create_flow.hpp"


int main() {
    //big_veg_derivs_ex();
    //calc_derivs_ex();

    constexpr output_settings o{};

    solve_flow<1000, 128, 128, 128>(o);



    return 0;
}
