//
// Created by jacob on 11/8/21.
//

#ifndef CODE_CALC_HPP
#define CODE_CALC_HPP

#include "vec3.hpp"
#include "big_vec.hpp"  //to removed later




//nx,ny,nz denote the derivative. e.g. <1,1,0> corresponds to d^2/dxdy
//ind corresponds to where the derivative is to take place
template<unsigned nx, unsigned ny, unsigned nz, typename T>
auto smart_deriv(const T& v, const unsigned ind) noexcept  {
    static_assert( (nx < 4) && (ny < 4) && (nz < 4), "Forth derivatives and higher are not supported");
    static_assert( (nx+ny+nz) < 4, "Mixed derivatives higher than third order are not supported" );
    //static_assert( !( (nx==1) && (ny==1) && (nz==1) ), "d^3/dxdydz is not currently supported" );
    static_assert( (nx+ny+nz) > 0, "Need at least 1 dimension for a derivative" );

    constexpr auto m = std::max( {nx, ny, nz} );
    constexpr auto s = nx + ny + nz;

    if constexpr ((nx==1) && (ny==1) && (nz==1)) {  //d^3/dxdydz
        const auto dx = v.dx(ind);
        //const auto dy = v.dy(ind);
        //const auto dz = v.dz(ind);

        const bool cx = v.can_move(ind, 1,0,0) && v.can_move(ind, -1,0,0);
        //const bool cy = v.can_move(0,1,0) && v.can_move(0,-1,0);
        //const bool cz = v.can_move(0,0,2) && v.can_move(0,0,-1);

        const bool bx = v.can_move(ind, -2,0,0);
        //const bool by = v.can_move(0,-2,0);
        //const bool bz = v.can_move(0,0,-2);

        //const bool fx = v.can_move(2,0,0);
        //const bool fy = v.can_move(0,2,0);
        //const bool fz = v.can_move(0,0,2);

        if (cx) {
            return 1/(2*dx) * ( smart_deriv<0,1,1>(v, v.get_move_ind(ind, 1,0,0)) - smart_deriv<0,1,1>(v, v.get_move_ind(ind, -1,0,0))  );
        } else if (bx) {
            return 1/(2*dx) * ( 3*smart_deriv<0,1,1>(v, ind) - 4*smart_deriv<0,1,1>(v, v.get_move_ind(ind, -1,0,0)) + smart_deriv<0,1,1>(v, v.get_move_ind(ind, -2,0,0))  );
        } else {
            return 1/(2*dx) * (-3*smart_deriv<0,1,1>(v, ind) + 4*smart_deriv<0,1,1>(v, v.get_move_ind(ind, 1,0,0))  - smart_deriv<0,1,1>(v, v.get_move_ind(ind,  2,0,0))  );
        }


    }

    if constexpr (s==1) {   //first order derivative
        constexpr unsigned axis = nx == 1 ? 0 : ny == 1 ? 1 : 2;

        constexpr vec3 move_vec1 = make_vec<axis, 1>();
        constexpr vec3 move_vec2 = make_vec<axis, -1>();
        constexpr vec3 move_vec3 = make_vec<axis, 2>();

        if (v.can_move(ind, move_vec1) && v.can_move(ind, move_vec2) ) {
            return central_difference_1st<axis>(v, ind);
        } else if (v.can_move(ind, move_vec3) ) {
            return forward_difference_1st<axis>(v, ind);
        } else {
            return backward_difference_1st<axis>(v, ind);
        }

    }
    else if constexpr(s==2) { //second order derivative
        if constexpr(m==2) {    //pure derivative
            constexpr unsigned axis = nx == 2 ? 0 : ny == 2 ? 1 : 2;

            constexpr vec3 move_vec1 = make_vec<axis, 1>();
            constexpr vec3 move_vec2 = make_vec<axis, -1>();
            constexpr vec3 move_vec3 = make_vec<axis, 3>();

            if (v.can_move(ind, move_vec1) && v.can_move(ind, move_vec2) ) {
                return central_difference_2nd<axis>(v, ind);
            } else if ( v.can_move(ind, move_vec3) ) {
                return forward_difference_2nd<axis>(v, ind);
            } else {
                return backward_difference_2nd<axis>(v, ind);
            }

        }
        else {    //mixed derivative
            constexpr unsigned axis1 = nx == 1 ? 0 : ny == 1 ? 1 : 2;
            constexpr unsigned axis2 = (ny == 1 && axis1!=1) ? 1 : 2;

            constexpr vec3 move_vec_c1 = make_vec<axis1, 1, axis2, 1>();
            constexpr vec3 move_vec_c2 = make_vec<axis1, -1, axis2, 1>();
            constexpr vec3 move_vec_c3 = make_vec<axis1, 1, axis2, -1>();
            constexpr vec3 move_vec_c4 = make_vec<axis1, -1, axis2, -1>();

            constexpr vec3 move_vec_f1_1 = make_vec<axis1, 1, axis2, 2>();
            //constexpr vec3 move_vec_f2_1 = make_vec<axis1, 1, axis2, 1>();
            //constexpr vec3 move_vec_f3_1 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_f4_1 = make_vec<axis1, -1, axis2, 2>();
            //constexpr vec3 move_vec_f5_1 = make_vec<axis1, -1, axis2, 1>();
            //constexpr vec3 move_vec_f6_1 = make_vec<axis1, -1, axis2, 0>();

            constexpr vec3 move_vec_f1_2 = make_vec<axis2, 1, axis1, 2>();
            constexpr vec3 move_vec_f2_2 = make_vec<axis2, 1, axis1, 1>();
            constexpr vec3 move_vec_f3_2 = make_vec<axis2, 1, axis1, 0>();
            constexpr vec3 move_vec_f4_2 = make_vec<axis2, -1, axis1, 2>();
            constexpr vec3 move_vec_f5_2 = make_vec<axis2, -1, axis1, 1>();
            constexpr vec3 move_vec_f6_2 = make_vec<axis2, -1, axis1, 0>();

            constexpr vec3 move_vec_b1_1 = make_vec<axis1, 1, axis2, -2>();
            //constexpr vec3 move_vec_b2_1 = make_vec<axis1, 1, axis2, -1>();
            //constexpr vec3 move_vec_b3_1 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_b4_1 = make_vec<axis1, -1, axis2, -2>();
            //constexpr vec3 move_vec_b5_1 = make_vec<axis1, -1, axis2, -1>();
            //constexpr vec3 move_vec_b6_1 = make_vec<axis1, -1, axis2, 0>();

            constexpr vec3 move_vec_b1_2 = make_vec<axis2, 1, axis1, -2>();
            constexpr vec3 move_vec_b2_2 = make_vec<axis2, 1, axis1, -1>();
            constexpr vec3 move_vec_b3_2 = move_vec_f3_2;//make_vec<axis2, 1, axis1, 0>();
            constexpr vec3 move_vec_b4_2 = make_vec<axis2, -1, axis1, -2>();
            constexpr vec3 move_vec_b5_2 = make_vec<axis2, -1, axis1, -1>();
            constexpr vec3 move_vec_b6_2 = move_vec_f6_2; //make_vec<axis2, -1, axis1, 0>();

            constexpr vec3 move_vec_fb1_1 = make_vec<axis1, 2, axis2, -2>();
            constexpr vec3 move_vec_fb2_1 = make_vec<axis1, 1, axis2, -2>();
            constexpr vec3 move_vec_fb3_1 = make_vec<axis1, 0, axis2, -2>();
            //constexpr vec3 move_vec_fb4_1 = make_vec<axis1, 2, axis2, -1>();
            //constexpr vec3 move_vec_fb5_1 = make_vec<axis1, 1, axis2, -1>();
            //constexpr vec3 move_vec_fb6_1 = make_vec<axis1, 0, axis2, -1>();
            //constexpr vec3 move_vec_fb7_1 = make_vec<axis1, 2, axis2, 0>();
            //constexpr vec3 move_vec_fb8_1 = make_vec<axis1, 1, axis2, 0>();

            constexpr vec3 move_vec_fb1_2 = make_vec<axis2, 2, axis1, -2>();
            //constexpr vec3 move_vec_fb2_2 = make_vec<axis2, 2, axis1, -1>();
            constexpr vec3 move_vec_fb3_2 = make_vec<axis2, 2, axis1, 0>();
            constexpr vec3 move_vec_fb4_2 = make_vec<axis2, 1, axis1, -2>();
            //constexpr vec3 move_vec_fb5_2 = make_vec<axis2, 1, axis1, -1>();
            //constexpr vec3 move_vec_fb6_2 = make_vec<axis2, 1, axis1, 0>();
            //constexpr vec3 move_vec_fb7_2 = make_vec<axis2, 0, axis1, -1>();
            constexpr vec3 move_vec_fb8_2 = make_vec<axis2, 0, axis1, -2>();

            constexpr vec3 move_vec_ff1 = make_vec<axis1, 2, axis2, 2>();
            constexpr vec3 move_vec_ff2 = make_vec<axis1, 1, axis2, 2>();
            constexpr vec3 move_vec_ff3 = make_vec<axis1, 0, axis2, 2>();
            //constexpr vec3 move_vec_ff4 = make_vec<axis1, 2, axis2, 1>();
            //constexpr vec3 move_vec_ff5 = make_vec<axis1, 1, axis2, 1>();
            //constexpr vec3 move_vec_ff6 = make_vec<axis1, 0, axis2, 1>();
            //constexpr vec3 move_vec_ff7 = make_vec<axis1, 2, axis2, 0>();
            //constexpr vec3 move_vec_ff8 = make_vec<axis1, 1, axis2, 0>();

            /*constexpr vec3 move_vec_bb1 = make_vec<axis1, -2, axis2, -2>();
            constexpr vec3 move_vec_bb2 = make_vec<axis1, -1, axis2, -2>();
            constexpr vec3 move_vec_bb3 = make_vec<axis1, 0, axis2, -2>();*/
            //constexpr vec3 move_vec_bb4 = make_vec<axis1, -2, axis2, -1>();
            //constexpr vec3 move_vec_bb5 = make_vec<axis1, -1, axis2, -1>();
            //constexpr vec3 move_vec_bb6 = make_vec<axis1, 0, axis2, -1>();
            //constexpr vec3 move_vec_bb7 = make_vec<axis1, -2, axis2, 0>();
            //constexpr vec3 move_vec_bb8 = make_vec<axis1, -1, axis2, 0>();



            if (v.can_move(ind, move_vec_c1) && v.can_move(ind, move_vec_c2)
            && v.can_move(ind, move_vec_c3) && v.can_move(ind, move_vec_c4)) {
                return central_difference_2nd_mixed<axis1, axis2>(v, ind);
            } else if ( v.can_move(ind, move_vec_f1_1) && v.can_move(ind, move_vec_f4_1) ) {
                return central_forward_difference_2nd_mixed<axis1, axis2>(v, ind);
            }else if ( v.can_move(ind, move_vec_f1_2) && v.can_move(ind, move_vec_f2_2) && v.can_move(ind, move_vec_f3_2)
             && v.can_move(ind, move_vec_f4_2) && v.can_move(ind, move_vec_f5_2) && v.can_move(ind, move_vec_f6_2)) {
                return central_forward_difference_2nd_mixed<axis2, axis1>(v, ind);
            } else if ( v.can_move(ind, move_vec_b1_1) && v.can_move(ind, move_vec_b4_1) ) {
                return central_backward_difference_2nd_mixed<axis1, axis2>(v, ind);
            }else if ( v.can_move(ind, move_vec_b1_2) && v.can_move(ind, move_vec_b2_2) && v.can_move(ind, move_vec_b3_2)
            && v.can_move(ind, move_vec_b4_2) && v.can_move(ind, move_vec_b5_2) && v.can_move(ind, move_vec_b6_2)) {
                return central_backward_difference_2nd_mixed<axis2, axis1>(v, ind);
            } else if ( v.can_move(ind, move_vec_fb1_1) && v.can_move(ind, move_vec_fb2_1)  && v.can_move(ind, move_vec_fb3_1) ) {
                return forward_backward_difference_2nd_mixed<axis1, axis2>(v, ind);
            } else if ( v.can_move(ind, move_vec_fb1_2) && v.can_move(ind, move_vec_fb3_2)
            && v.can_move(ind, move_vec_fb4_2)  && v.can_move(ind, move_vec_fb8_2) ) {
                return forward_backward_difference_2nd_mixed<axis2, axis1>(v, ind);
            } else if ( v.can_move(ind, move_vec_ff1) && v.can_move(ind, move_vec_ff2)  && v.can_move(ind, move_vec_ff3) ) {
                return forward_difference_2nd_mixed<axis1, axis2>(v, ind);
            } else {//if ( v.can_move(ind, move_vec_bb1) && v.can_move(ind, move_vec_bb2)  && v.can_move(ind, move_vec_bb3) ) {
                return backward_difference_2nd_mixed<axis1, axis2>(v, ind);
            }


        }

    }
    else if constexpr(s==3) { //third order derivative
        if constexpr(m==3) {    //pure derivative
            constexpr unsigned axis = nx == 3 ? 0 : ny == 3 ? 1 : 2;

            constexpr vec3 move_vec_c1 = make_vec<axis, 2>();
            constexpr vec3 move_vec_c2 = make_vec<axis, -2>();
            constexpr vec3 move_vec_f = make_vec<axis, 4>();

            if (v.can_move(ind, move_vec_c1) && v.can_move(ind, move_vec_c2)) {
                return central_difference_3rd<axis>(v, ind);
            } else if (v.can_move(ind, move_vec_f)) {
                return forward_difference_3rd<axis>(v, ind);
            } else {
                return backward_difference_3rd<axis>(v, ind);
            }

        } else {    //mixed derivative
            constexpr unsigned axis1 = nx == 1 ? 0 : ny == 1 ? 1 : 2;
            constexpr unsigned axis2 = nx == 2 ? 0 : ny == 2 ? 1 : 2;

            constexpr vec3 move_vec_cc1 = make_vec<axis1, 1, axis2, 1>();
            constexpr vec3 move_vec_cc2 = make_vec<axis1, 0, axis2, 1>();
            constexpr vec3 move_vec_cc3 = make_vec<axis1, -1, axis2, 1>();
            constexpr vec3 move_vec_cc4 = make_vec<axis1, 1, axis2, -1>();
            constexpr vec3 move_vec_cc5 = make_vec<axis1, 0, axis2, -1>();
            constexpr vec3 move_vec_cc6 = make_vec<axis1, -1, axis2, -1>();

            constexpr vec3 move_vec_ff1 = make_vec<axis1, 2, axis2, 3>();
            constexpr vec3 move_vec_ff2 = make_vec<axis1, 2, axis2, 2>();
            constexpr vec3 move_vec_ff3 = make_vec<axis1, 2, axis2, 1>();
            constexpr vec3 move_vec_ff4 = make_vec<axis1, 2, axis2, 0>();
            constexpr vec3 move_vec_ff5 = make_vec<axis1, 1, axis2, 3>();
            constexpr vec3 move_vec_ff6 = make_vec<axis1, 1, axis2, 2>();
            constexpr vec3 move_vec_ff7 = make_vec<axis1, 1, axis2, 1>();
            constexpr vec3 move_vec_ff8 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_ff9 = make_vec<axis1, 0, axis2, 3>();
            constexpr vec3 move_vec_ff10 = make_vec<axis1, 0, axis2, 2>();
            constexpr vec3 move_vec_ff11 = make_vec<axis1, 0, axis2, 1>();

            constexpr vec3 move_vec_bb1 = make_vec<axis1, -2, axis2, -3>();
            constexpr vec3 move_vec_bb2 = make_vec<axis1, -2, axis2, -2>();
            constexpr vec3 move_vec_bb3 = make_vec<axis1, -2, axis2, -1>();
            constexpr vec3 move_vec_bb4 = make_vec<axis1, -2, axis2, -0>();
            constexpr vec3 move_vec_bb5 = make_vec<axis1, -1, axis2, -3>();
            constexpr vec3 move_vec_bb6 = make_vec<axis1, -1, axis2, -2>();
            constexpr vec3 move_vec_bb7 = make_vec<axis1, -1, axis2, -1>();
            constexpr vec3 move_vec_bb8 = make_vec<axis1, -1, axis2, 0>();
            constexpr vec3 move_vec_bb9 = make_vec<axis1, 0, axis2, -3>();
            constexpr vec3 move_vec_bb10 = make_vec<axis1, 0, axis2, -2>();
            constexpr vec3 move_vec_bb11 = make_vec<axis1, 0, axis2, -1>();

            constexpr vec3 move_vec_cf1 = make_vec<axis1, 2, axis2, 1>();
            constexpr vec3 move_vec_cf2 = make_vec<axis1, 2, axis2, 0>();
            constexpr vec3 move_vec_cf3 = make_vec<axis1, 2, axis2, -1>();
            constexpr vec3 move_vec_cf4 = make_vec<axis1, 1, axis2, 1>();
            constexpr vec3 move_vec_cf5 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_cf6 = make_vec<axis1, 1, axis2, -1>();
            constexpr vec3 move_vec_cf7 = make_vec<axis1, 0, axis2, 1>();
            constexpr vec3 move_vec_cf8 = make_vec<axis1, 0, axis2, -1>();

            constexpr vec3 move_vec_cb1 = make_vec<axis1, -2, axis2, 1>();
            constexpr vec3 move_vec_cb2 = make_vec<axis1, -2, axis2, 0>();
            constexpr vec3 move_vec_cb3 = make_vec<axis1, -2, axis2, -1>();
            constexpr vec3 move_vec_cb4 = make_vec<axis1, -1, axis2, 1>();
            constexpr vec3 move_vec_cb5 = make_vec<axis1, -1, axis2, 0>();
            constexpr vec3 move_vec_cb6 = make_vec<axis1, -1, axis2, -1>();
            constexpr vec3 move_vec_cb7 = make_vec<axis1, 0, axis2, 1>();
            constexpr vec3 move_vec_cb8 = make_vec<axis1, 0, axis2, -1>();

            constexpr vec3 move_vec_fc1 = make_vec<axis1, 1, axis2, 3>();
            constexpr vec3 move_vec_fc2 = make_vec<axis1, 1, axis2, 2>();
            constexpr vec3 move_vec_fc3 = make_vec<axis1, 1, axis2, 1>();
            constexpr vec3 move_vec_fc4 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_fc5 = make_vec<axis1, -1, axis2, 3>();
            constexpr vec3 move_vec_fc6 = make_vec<axis1, -1, axis2, 2>();
            constexpr vec3 move_vec_fc7 = make_vec<axis1, -1, axis2, 1>();
            constexpr vec3 move_vec_fc8 = make_vec<axis1, -1, axis2, 0>();

            constexpr vec3 move_vec_fb1 = make_vec<axis1, -2, axis2, 3>();
            constexpr vec3 move_vec_fb2 = make_vec<axis1, -2, axis2, 2>();
            constexpr vec3 move_vec_fb3 = make_vec<axis1, -2, axis2, 1>();
            constexpr vec3 move_vec_fb4 = make_vec<axis1, -2, axis2, 0>();
            constexpr vec3 move_vec_fb5 = make_vec<axis1, -1, axis2, 3>();
            constexpr vec3 move_vec_fb6 = make_vec<axis1, -1, axis2, 2>();
            constexpr vec3 move_vec_fb7 = make_vec<axis1, -1, axis2, 1>();
            constexpr vec3 move_vec_fb8 = make_vec<axis1, -1, axis2, 0>();
            constexpr vec3 move_vec_fb9 = make_vec<axis1, 0, axis2, 3>();
            constexpr vec3 move_vec_fb10 = make_vec<axis1, 0, axis2, 2>();
            constexpr vec3 move_vec_fb11 = make_vec<axis1, 0, axis2, 1>();

            constexpr vec3 move_vec_bc1 = make_vec<axis1, 1, axis2, -3>();
            constexpr vec3 move_vec_bc2 = make_vec<axis1, 1, axis2, -2>();
            constexpr vec3 move_vec_bc3 = make_vec<axis1, 1, axis2, -1>();
            constexpr vec3 move_vec_bc4 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_bc5 = make_vec<axis1, -1, axis2, -3>();
            constexpr vec3 move_vec_bc6 = make_vec<axis1, -1, axis2, -2>();
            constexpr vec3 move_vec_bc7 = make_vec<axis1, -1, axis2, -1>();
            constexpr vec3 move_vec_bc8 = make_vec<axis1, -1, axis2, 0>();

            constexpr vec3 move_vec_bf1 = make_vec<axis1, 2, axis2, -3>();
            constexpr vec3 move_vec_bf2 = make_vec<axis1, 2, axis2, -2>();
            constexpr vec3 move_vec_bf3 = make_vec<axis1, 2, axis2, -1>();
            constexpr vec3 move_vec_bf4 = make_vec<axis1, 2, axis2, 0>();
            constexpr vec3 move_vec_bf5 = make_vec<axis1, 1, axis2, -3>();
            constexpr vec3 move_vec_bf6 = make_vec<axis1, 1, axis2, -2>();
            constexpr vec3 move_vec_bf7 = make_vec<axis1, 1, axis2, -1>();
            constexpr vec3 move_vec_bf8 = make_vec<axis1, 1, axis2, 0>();
            constexpr vec3 move_vec_bf9 = make_vec<axis1, 0, axis2, -3>();
            constexpr vec3 move_vec_bf10 = make_vec<axis1, 0, axis2, -2>();
            constexpr vec3 move_vec_bf11 = make_vec<axis1, 0, axis2, -1>();

            if (v.can_move(ind, move_vec_cc1) && v.can_move(ind, move_vec_cc2) && v.can_move(ind, move_vec_cc3) &&
            v.can_move(ind, move_vec_cc4) && v.can_move(ind, move_vec_cc5) && v.can_move(ind, move_vec_cc6) ) {
                return central_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_cf1) && v.can_move(ind, move_vec_cf2) && v.can_move(ind, move_vec_cf3) &&
            v.can_move(ind, move_vec_cf4) && v.can_move(ind, move_vec_cf5) && v.can_move(ind, move_vec_cf6) &&
            v.can_move(ind, move_vec_cf7) && v.can_move(ind, move_vec_cf8) ) {
                return central_forward_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_cb1) && v.can_move(ind, move_vec_cb2) && v.can_move(ind, move_vec_cb3) &&
            v.can_move(ind, move_vec_cb4) && v.can_move(ind, move_vec_cb5) && v.can_move(ind, move_vec_cb6) &&
            v.can_move(ind, move_vec_cb7) && v.can_move(ind, move_vec_cb8) ) {
                return central_backward_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_fc1) && v.can_move(ind, move_vec_fc2) && v.can_move(ind, move_vec_fc3) &&
            v.can_move(ind, move_vec_fc4) && v.can_move(ind, move_vec_fc5) && v.can_move(ind, move_vec_fc6) &&
            v.can_move(ind, move_vec_fc7) && v.can_move(ind, move_vec_fc8)  ) {
                return forward_central_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_fb1) && v.can_move(ind, move_vec_fb2) && v.can_move(ind, move_vec_fb3) &&
            v.can_move(ind, move_vec_fb4) && v.can_move(ind, move_vec_fb5) && v.can_move(ind, move_vec_fb6) &&
            v.can_move(ind, move_vec_fb7) && v.can_move(ind, move_vec_fb8) && v.can_move(ind, move_vec_fb9) &&
            v.can_move(ind, move_vec_fb10) && v.can_move(ind, move_vec_fb11) ) {
                return forward_backward_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_bc1) && v.can_move(ind, move_vec_bc2) && v.can_move(ind, move_vec_bc3) &&
            v.can_move(ind, move_vec_bc4) && v.can_move(ind, move_vec_bc5) && v.can_move(ind, move_vec_bc6) &&
            v.can_move(ind, move_vec_bc7) && v.can_move(ind, move_vec_bc8) ) {
                return backward_central_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_ff1) && v.can_move(ind, move_vec_ff2) && v.can_move(ind, move_vec_ff3) &&
              v.can_move(ind, move_vec_ff4) && v.can_move(ind, move_vec_ff5) && v.can_move(ind, move_vec_ff6) &&
              v.can_move(ind, move_vec_ff7) && v.can_move(ind, move_vec_ff8) && v.can_move(ind, move_vec_ff9) &&
              v.can_move(ind, move_vec_ff10) && v.can_move(ind, move_vec_ff11)) {
                return forward_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else if (v.can_move(ind, move_vec_bb1) && v.can_move(ind, move_vec_bb2) && v.can_move(ind, move_vec_bb3) &&
            v.can_move(ind, move_vec_bb4) && v.can_move(ind, move_vec_bb5) && v.can_move(ind, move_vec_bb6) &&
            v.can_move(ind, move_vec_bb7) && v.can_move(ind, move_vec_bb8) && v.can_move(ind, move_vec_bb9) &&
            v.can_move(ind, move_vec_bb10) && v.can_move(ind, move_vec_bb11) ) {
                return backward_difference_3rd_mixed<axis2, axis1>(v, ind);
            } else { //if (v.can_move(ind, move_vec_bf1) && v.can_move(ind, move_vec_bf2) && v.can_move(ind, move_vec_bf3) &&
            //v.can_move(ind, move_vec_bf4) && v.can_move(ind, move_vec_bf5) && v.can_move(ind, move_vec_bf6) &&
            //v.can_move(ind, move_vec_bf7) && v.can_move(ind, move_vec_bf8) && v.can_move(ind, move_vec_bf9) &&
            //v.can_move(ind, move_vec_bf10) && v.can_move(ind, move_vec_bf11) ) {
                return backward_forward_difference_3rd_mixed<axis2, axis1>(v, ind);
            }

        }
    }
}


//H(v_{ind})
vec3 advection(const big_vec_v& v, const unsigned ind) noexcept {
    const auto s_v = v(ind);  //small v
    return s_v.x() * smart_deriv<1,0,0>(v, ind) + s_v.y() * smart_deriv<0,1,0>(v, ind) + s_v.z() * smart_deriv<0,0,1>(v, ind) ;
}

//L(v_{ind})
template <typename T>
auto laplacian(const T& v, const unsigned ind) noexcept {
    return smart_deriv<2,0,0>(v, ind) + smart_deriv<0,2,0>(v, ind) +smart_deriv<0,0,2>(v, ind);
}

//TODO : Test
vec3 gradient(const big_vec_d& v, const unsigned ind) noexcept {
    return vec3(smart_deriv<1,0,0>(v, ind), smart_deriv<0,1,0>(v, ind), smart_deriv<0,0,1>(v, ind));
}

//D(v_{ind})
double divergence(const big_vec_v& v, const unsigned ind) noexcept  {
    return smart_deriv<1,0,0>(v.xv, ind) + smart_deriv<0,1,0>(v.yv, ind) + smart_deriv<0,0,1>(v.zv, ind);
}

//D(H(v_{ind}))
double divergence_advection(const big_vec_v& v, const unsigned ind) noexcept {
    const auto s_v = v(ind);  //small vec

    const auto vx_x = smart_deriv<1,0,0>(v.xv,ind);
    const auto vy_x = smart_deriv<0,1,0>(v.xv,ind);
    const auto vz_x = smart_deriv<0,0,1>(v.xv,ind);

    const auto vx_y = smart_deriv<1,0,0>(v.yv,ind);
    const auto vy_y = smart_deriv<0,1,0>(v.yv,ind);
    const auto vz_y = smart_deriv<0,0,1>(v.yv,ind);

    const auto vx_z = smart_deriv<1,0,0>(v.zv,ind);
    const auto vy_z = smart_deriv<0,1,0>(v.zv,ind);
    const auto vz_z = smart_deriv<0,0,1>(v.zv,ind);


    const auto vxx_x = smart_deriv<2,0,0>(v.xv,ind);
    const auto vyy_y = smart_deriv<0,2,0>(v.yv,ind);
    const auto vzz_z = smart_deriv<0,0,2>(v.zv,ind);


    const auto vxy_x = smart_deriv<1,1,0>(v.xv,ind);

    const auto vxz_x = smart_deriv<1,0,1>(v.xv,ind);

    const auto vxy_y = smart_deriv<1,1,0>(v.yv,ind);
    const auto vyz_y = smart_deriv<0,1,1>(v.yv,ind);

    const auto vxz_z = smart_deriv<1,0,1>(v.zv,ind);
    const auto vyz_z = smart_deriv<0,1,1>(v.zv,ind);

    const auto term1 = vx_x*vx_x + vy_y*vy_y + vz_z*vz_z;
    const auto term2 = 2*vx_y*vy_x + 2*vx_z*vz_x + 2*vy_z*vz_y;
    const auto term3 = s_v.x() * ( vxx_x + vxy_y + vxz_z );
    const auto term4 = s_v.y() * ( vxy_x + vyy_y + vyz_z );
    const auto term5 = s_v.z() * ( vxz_x + vyz_y + vzz_z );


    return term1 + term2 + term3 + term4 + term5;
}

//D(L(v_{ind}))
double divergence_laplacian(const big_vec_v& v, const unsigned ind) noexcept {
    const auto vxxx = smart_deriv<3,0,0>(v.xv,ind);
    const auto vxyy = smart_deriv<1,2,0>(v.xv,ind);
    const auto vxzz = smart_deriv<1,0,2>(v.xv,ind);

    const auto vyxx = smart_deriv<2,1,0>(v.yv,ind);
    const auto vyyy = smart_deriv<0,3,0>(v.yv,ind);
    const auto vyzz = smart_deriv<0,1,2>(v.yv,ind);

    const auto vzxx = smart_deriv<2,0,1>(v.zv,ind);
    const auto vzyy = smart_deriv<0,2,1>(v.zv,ind);
    const auto vzzz = smart_deriv<0,0,3>(v.zv,ind);

    return vxxx+ vxyy + vxzz + vyxx + vyyy + vyzz + vzxx + vzyy + vzzz;
}




#endif //CODE_CALC_HPP
