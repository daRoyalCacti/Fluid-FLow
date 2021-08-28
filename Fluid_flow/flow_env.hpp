//
// Created by jacob on 21/8/21.
//

#ifndef CODE_FLOW_ENV_HPP
#define CODE_FLOW_ENV_HPP

#include "../Rigid_body/mesh.hpp"
#include "../MyMath/big_vec.hpp"
#include "../Rigid_body/triangle.hpp"



template <unsigned N, unsigned M, unsigned P>
void v_IC(big_vec<N,M,P, vec3> &v) noexcept {
    //just leave it as 0
    for (unsigned i = 0; i <= N; i++) {
        for (unsigned j = 0; j <= M; j++) {
            for (unsigned k = 0; k <= P; k++) {
                v.add_elm(i, j, k, 0, 0, 0);
            }
        }

    }
}


#endif //CODE_FLOW_ENV_HPP
