//
// Created by jacob on 21/8/21.
//

#ifndef CODE_FLOW_ENV_HPP
#define CODE_FLOW_ENV_HPP

#include "../Rigid_body/mesh.hpp"
#include "../MyMath/big_vec.hpp"
#include "../Rigid_body/triangle.hpp"



void v_IC(big_vec_v &v) noexcept {
    //just leave it as 0
    for (unsigned i = 0; i < v.size(); i++) {
        v.add_elm(i, 0,0,0);
    }

}


#endif //CODE_FLOW_ENV_HPP
