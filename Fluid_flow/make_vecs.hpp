//
// Created by jacob on 12/8/21.
//

#ifndef CODE_MAKE_VECS_HPP
#define CODE_MAKE_VECS_HPP

#include "../MyMath/calc.hpp"

//TODO : test to see if gives right output
void make_b(big_vec_v &b, const double Re, const double dt, const big_vec_v &v_n, const big_vec_v &v_n1, const big_vec_d &p, const boundary_conditions &bc) noexcept {
#pragma omp parallel for
    for (unsigned i = 0; i < p.size(); i++) {
        if (v_n.is_boundary(i)) {
            if (v_n.g->off_walls(i)) {
                b.add_elm(i,  bc.v_points.get_vel(i) - v_n(i)); //need to minus v_n(i) because solving for pressure correction
            } else {
                b.add_elm(i, 0,0,0);
            }
        } else {
#ifdef VC
            b.add_elm(i,
                      -3 / 2 * advection(v_n, i) + 1 / 2 * advection(v_n1, i)
                      - gradient(p, i) + 1 / Re * laplacian(v_n, i) );
#else
            b.add_elm(i,
                      -3 / 2 * advection(v_n, i) + 1 / (2 * Re) * laplacian(v_n, i) +
                      v_n(i) / dt -
                      gradient(p, i) + 1 / 2 * advection(v_n1, i));
#endif
        }

    }


}

//TODO : test to see if gives right output
//for the first timestep
void make_b_first(big_vec_v &b, const double Re, const double dt, const big_vec_v &v_n, const big_vec_d &p, const boundary_conditions &bc) noexcept {
#pragma omp parallel for
    for (unsigned i = 0; i < p.size(); i++) {
        if (v_n.is_boundary(i)) {
            if (v_n.g->off_walls(i)) {
                b.add_elm(i,  bc.v_points.get_vel(i) - v_n(i)); //need to minus v_n(i) because solving for pressure correction
            } else {
                b.add_elm(i, 0,0,0);
            }
        } else {
#ifdef VC
            b.add_elm(i,
                      -advection(v_n, i) + 1 / (Re) * laplacian(v_n, i) - gradient(p, i));
#else
            b.add_elm(i,
                      -advection(v_n, i) + 1 / (2 * Re) * laplacian(v_n, i) + v_n(i) / dt -
                      gradient(p, i));
#endif
        }
    }

}

//TODO : test to see if gives right output
void make_s(big_vec_d &s, const double Re, const double dt, const big_vec_v &v_n,
        const big_vec_v &v_n1, const big_vec_d &p) noexcept {
#pragma omp parallel for
    for (unsigned i = 0; i < p.size(); i++) {
        if (p.is_boundary(i)) {
            if (p.g->off_walls(i)) {
                s(i) = 0;
            } else {
                s(i) = 0-p(i);  //need to minus p(i) because solving for pressure correction
            }
            //s(i) = p(i);
        } else {

            s(i) = divergence(v_n, i) / dt - 3/2 * divergence_advection(v_n, i) + 1/2 *
                    divergence_advection(v_n1, i) + 3/(2*Re) *
                    divergence_laplacian(v_n,i) - 1/(2*Re) * divergence_laplacian(v_n1, i) - laplacian(p, i);

        }



    }
}


//TODO : test to see if gives right output
void make_s_first(big_vec_d &s, const double Re, const double dt, const big_vec_v &v_n, const big_vec_d &p) noexcept {
    #pragma omp parallel for
    for (unsigned i = 0; i < p.size(); i++) {
        if (p.is_boundary(i)) {
            if (p.g->off_walls(i)) {
                s(i) = 0;
            } else {
                s(i) = 0-p(i);  //need to minus p(i) because solving for pressure correction
            }
            //s(i) = p(i);
        } else {
            s(i) = divergence(v_n, i) / dt - divergence_advection(v_n, i) + 1 / Re *
                    divergence_laplacian(v_n, i) - laplacian(p, i);
        }

    }
}

#endif //CODE_MAKE_VECS_HPP
