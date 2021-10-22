//
// Created by jacob on 11/8/21.
//

#ifndef CODE_CALC_HPP
#define CODE_CALC_HPP

#include "big_vec.hpp"




//nx,ny,nz denote the derivative. e.g. <1,1,0> corresponds to d^2/dxdy
//ind corresponds to where the derivative is to take place
template<unsigned nx, unsigned ny, unsigned nz, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
        inline auto smart_deriv_old(const T& v, const unsigned ind) noexcept  {
    static_assert( (nx < 4) && (ny < 4) && (nz < 4), "Forth derivatives and higher are not supported");
    static_assert( (nx+ny+nz) < 4, "Mixed derivatives higher than third order are not supported" );
    static_assert( (nx+ny+nz) > 0, "Need at least 1 dimension for a derivative" );

    constexpr auto m = std::max( {nx, ny, nz} );
    constexpr auto s = nx + ny + nz;

    if constexpr ((nx==1) && (ny==1) && (nz==1)) {  //d^3/dxdydz

        if (v.can_move(ind, 1,0,0) && v.can_move(ind, -1,0,0)) {    //central
            return 1/(2*v.dx(ind)) * ( smart_deriv_old<0,1,1>(v, v.get_move_ind(ind, 1,0,0)) - smart_deriv_old<0,1,1>(v, v.get_move_ind(ind, -1,0,0))  );
        } else if (v.can_move(ind, -2,0,0)) {   //backward
            return 1/(2*v.dx(ind)) * ( 3*smart_deriv_old<0,1,1>(v, ind) - 4*smart_deriv_old<0,1,1>(v, v.get_move_ind(ind, -1,0,0)) + smart_deriv_old<0,1,1>(v, v.get_move_ind(ind, -2,0,0))  );
        } else {
            return 1/(2*v.dx(ind)) * (-3*smart_deriv_old<0,1,1>(v, ind) + 4*smart_deriv_old<0,1,1>(v, v.get_move_ind(ind, 1,0,0))  - smart_deriv_old<0,1,1>(v, v.get_move_ind(ind,  2,0,0))  );
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
            constexpr vec3 move_vec_f4_1 = make_vec<axis1, -1, axis2, 2>();

            constexpr vec3 move_vec_f1_2 = make_vec<axis2, 1, axis1, 2>();
            constexpr vec3 move_vec_f2_2 = make_vec<axis2, 1, axis1, 1>();
            constexpr vec3 move_vec_f3_2 = make_vec<axis2, 1, axis1, 0>();
            constexpr vec3 move_vec_f4_2 = make_vec<axis2, -1, axis1, 2>();
            constexpr vec3 move_vec_f5_2 = make_vec<axis2, -1, axis1, 1>();
            constexpr vec3 move_vec_f6_2 = make_vec<axis2, -1, axis1, 0>();

            constexpr vec3 move_vec_b1_1 = make_vec<axis1, 1, axis2, -2>();
            constexpr vec3 move_vec_b4_1 = make_vec<axis1, -1, axis2, -2>();

            constexpr vec3 move_vec_b1_2 = make_vec<axis2, 1, axis1, -2>();
            constexpr vec3 move_vec_b2_2 = make_vec<axis2, 1, axis1, -1>();
            constexpr vec3 move_vec_b3_2 = move_vec_f3_2;
            constexpr vec3 move_vec_b4_2 = make_vec<axis2, -1, axis1, -2>();
            constexpr vec3 move_vec_b5_2 = make_vec<axis2, -1, axis1, -1>();
            constexpr vec3 move_vec_b6_2 = move_vec_f6_2;

            constexpr vec3 move_vec_fb1_1 = make_vec<axis1, 2, axis2, -2>();
            constexpr vec3 move_vec_fb2_1 = make_vec<axis1, 1, axis2, -2>();
            constexpr vec3 move_vec_fb3_1 = make_vec<axis1, 0, axis2, -2>();

            constexpr vec3 move_vec_fb1_2 = make_vec<axis2, 2, axis1, -2>();
            constexpr vec3 move_vec_fb3_2 = make_vec<axis2, 2, axis1, 0>();
            constexpr vec3 move_vec_fb4_2 = make_vec<axis2, 1, axis1, -2>();
            constexpr vec3 move_vec_fb8_2 = make_vec<axis2, 0, axis1, -2>();

            constexpr vec3 move_vec_ff1 = make_vec<axis1, 2, axis2, 2>();
            constexpr vec3 move_vec_ff2 = make_vec<axis1, 1, axis2, 2>();
            constexpr vec3 move_vec_ff3 = make_vec<axis1, 0, axis2, 2>();



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
            } else {
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
            } else {
                return backward_forward_difference_3rd_mixed<axis2, axis1>(v, ind);
            }

        }
    }
}


template<unsigned nx, unsigned ny, unsigned nz, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
        inline auto smart_deriv(const axes& a, const T& v, const unsigned ind) noexcept  {
    static_assert( (nx < 4) && (ny < 4) && (nz < 4), "Forth derivatives and higher are not supported");
    static_assert( (nx+ny+nz) < 4, "Mixed derivatives higher than third order are not supported" );
    static_assert( (nx+ny+nz) > 0, "Need at least 1 dimension for a derivative" );

    if constexpr(nx + ny + nz == 1) {   //first order
        constexpr unsigned n = (nx==1) ? 0 : (ny==1) ? 1 : 2;
        const vec3& x = a.get<n>();
        return x.x()*smart_deriv_old<1,0,0>(v,ind) + x.y()*smart_deriv_old<0,1,0>(v,ind) + x.z()*smart_deriv_old<0,0,1>(v,ind);
    } else if constexpr(nx + ny + nz ==2) {  //second order
        constexpr unsigned highest_order = std::max({nx, ny, nz});
        vec3 x, y;
        if (highest_order == 1) {
            constexpr unsigned axis1 = (nx==1) ? 0 : (ny==1) ? 1 : 2;
            constexpr unsigned axis2 = (nx==1 && axis1!=0) ? 0 : (ny==1 && axis1!=1) ? 1 : 2;
            x = a.get<axis1>();
            y = a.get<axis2>();
        } else {
            constexpr unsigned axis = (nx==2) ? 0 : (ny==2) ? 1 : 2;
            x = a.get<axis>();
            y = a.get<axis>();
        }

        return x.x()*y.x()*smart_deriv_old<2,0,0>(v,ind) + x.y()*y.y()* smart_deriv_old<0,2,0>(v,ind) + x.z()*y.z()*
        smart_deriv_old<0,0,2>(v,ind) + (x.y()*y.x() + x.x()*y.y() )*smart_deriv_old<1,1,0>(v,ind) + (x.z()*y.x() + x.x()*y.z())*
        smart_deriv_old<1,0,1>(v,ind) + (x.z()*y.y() + x.y()*y.z())* smart_deriv_old<0,1,1>(v,ind);

    } else {    //third order
        static_assert( std::max({nx,ny,nz}) != 1, "d^3/dxdydz is not supported" );

        vec3 x,y;
        if constexpr(std::max({nx,ny,nz})==3) {
            constexpr unsigned axis = (nx==3) ? 0 : (ny==3) ? 1 : 2;
            x = a.get<axis>();
            y = a.get<axis>();
        } else {
            constexpr unsigned axis1 = (nx==2) ? 0 : (ny==2) ? 1 : 2;
            constexpr unsigned axis2 = (nx==1) ? 0 : (ny==1) ? 1 : 2;
            x = a.get<axis1>();
            y = a.get<axis2>();
        }

        return x.x()*x.x()*y.x()*smart_deriv_old<3,0,0>(v,ind) + (x.y()*x.y()*y.x() + 2*x.x()*x.y()*y.y())*smart_deriv_old<1,2,0>(v,ind) +
        (x.z()*x.z()*y.x() + 2*x.x()*x.z()*y.z())*smart_deriv_old<1,0,2>(v,ind) + (2*x.x()*x.y()*y.x() + x.x()*x.x()*y.y())*smart_deriv_old<2,1,0>(v,ind) +
        (2*x.x()*x.z()*y.x() + x.x()*x.x()*y.z())* smart_deriv_old<2,0,1>(v, ind) +
        (2*x.y()*x.z()*y.x() + 2*x.x()*x.z()*y.y() + 2*x.x()*x.y()*y.z())* smart_deriv_old<1,1,1>(v, ind) +
        x.y()*x.y()*y.y()*smart_deriv_old<0,3,0>(v,ind) + x.z()*x.z()*y.z()*smart_deriv_old<0,0,3>(v, ind) +
        (x.z()*x.z()*y.y() + 2*x.y()*x.z()*y.z())* smart_deriv_old<0,1,2>(v, ind) +
        (2*x.y()*x.z()*y.y() + x.y()*x.y()*y.z())* smart_deriv_old<0,2,1>(v, ind);

    }

}


//H(v_{ind})
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline vec3 advection_old(const big_vec_v& v, const unsigned ind) noexcept {
    const auto s_v = v(ind);  //small v
    return s_v.x() * smart_deriv_old<1,0,0>(v, ind) + s_v.y() * smart_deriv_old<0,1,0>(v, ind) + s_v.z() * smart_deriv_old<0,0,1>(v, ind) ;
}

#ifdef ALWAYS_INLINE_DERIVS
_attribute__((always_inline))
#endif
inline vec3 advection(const axes& a, const big_vec_v& v, const unsigned ind) noexcept {
    const auto s_v = v(ind);  //small v
    return s_v.x() * smart_deriv<1,0,0>(a, v, ind) + s_v.y() * smart_deriv<0,1,0>(a, v, ind) + s_v.z() * smart_deriv<0,0,1>(a, v, ind) ;
}

//L(v_{ind})
template <typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline auto laplacian_old(const T& v, const unsigned ind) noexcept {
    return smart_deriv_old<2,0,0>(v, ind) + smart_deriv_old<0,2,0>(v, ind) +smart_deriv_old<0,0,2>(v, ind);
}

template <typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline auto laplacian(const axes& a, const T& v, const unsigned ind) noexcept {
    return laplacian_old(v,ind);
}



#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline vec3 gradient_old(const big_vec_d& v, const unsigned ind) noexcept {
    return {smart_deriv_old<1,0,0>(v, ind), smart_deriv_old<0,1,0>(v, ind), smart_deriv_old<0,0,1>(v, ind)};
}

#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline vec3 gradient(const axes& a, const big_vec_d& v, const unsigned ind) noexcept {
    return {smart_deriv<1,0,0>(a, v, ind), smart_deriv<0,1,0>(a, v, ind), smart_deriv<0,0,1>(a, v, ind)};
}

//D(v_{ind})
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline double divergence_old(const big_vec_v& v, const unsigned ind) noexcept  {
    return smart_deriv_old<1,0,0>(v.xv, ind) + smart_deriv_old<0,1,0>(v.yv, ind) + smart_deriv_old<0,0,1>(v.zv, ind);
}

#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline double divergence(const axes& a, const big_vec_v& v, const unsigned ind) noexcept  {
    return divergence_old(v, ind);
    //return smart_deriv<1,0,0>(a, v.xv, ind) + smart_deriv<0,1,0>(a, v.yv, ind) + smart_deriv<0,0,1>(a, v.zv, ind);
}




//D(H(v_{ind}))
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline double divergence_advection_old(const big_vec_v& v, const unsigned ind) noexcept {
    const auto s_v = v(ind);  //small vec

    const auto vx_x = smart_deriv_old<1,0,0>(v.xv,ind);
    const auto vy_x = smart_deriv_old<0,1,0>(v.xv,ind);
    const auto vz_x = smart_deriv_old<0,0,1>(v.xv,ind);

    const auto vx_y = smart_deriv_old<1,0,0>(v.yv,ind);
    const auto vy_y = smart_deriv_old<0,1,0>(v.yv,ind);
    const auto vz_y = smart_deriv_old<0,0,1>(v.yv,ind);

    const auto vx_z = smart_deriv_old<1,0,0>(v.zv,ind);
    const auto vy_z = smart_deriv_old<0,1,0>(v.zv,ind);
    const auto vz_z = smart_deriv_old<0,0,1>(v.zv,ind);


    const auto vxx_x = smart_deriv_old<2,0,0>(v.xv,ind);
    const auto vyy_y = smart_deriv_old<0,2,0>(v.yv,ind);
    const auto vzz_z = smart_deriv_old<0,0,2>(v.zv,ind);


    const auto vxy_x = smart_deriv_old<1,1,0>(v.xv,ind);

    const auto vxz_x = smart_deriv_old<1,0,1>(v.xv,ind);

    const auto vxy_y = smart_deriv_old<1,1,0>(v.yv,ind);
    const auto vyz_y = smart_deriv_old<0,1,1>(v.yv,ind);

    const auto vxz_z = smart_deriv_old<1,0,1>(v.zv,ind);
    const auto vyz_z = smart_deriv_old<0,1,1>(v.zv,ind);

    const auto term1 = vx_x*vx_x + vy_y*vy_y + vz_z*vz_z;
    const auto term2 = 2*vx_y*vy_x + 2*vx_z*vz_x + 2*vy_z*vz_y;
    const auto term3 = s_v.x() * ( vxx_x + vxy_y + vxz_z );
    const auto term4 = s_v.y() * ( vxy_x + vyy_y + vyz_z );
    const auto term5 = s_v.z() * ( vxz_x + vyz_y + vzz_z );


    return term1 + term2 + term3 + term4 + term5;
}


#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline double divergence_advection(const axes& a, const big_vec_v& v, const unsigned ind) noexcept {
    const auto s_v = v(ind);  //small vec

    const auto vx_x = smart_deriv<1,0,0>(a, v.xv,ind);
    const auto vy_x = smart_deriv<0,1,0>(a, v.xv,ind);
    const auto vz_x = smart_deriv<0,0,1>(a, v.xv,ind);

    const auto vx_y = smart_deriv<1,0,0>(a, v.yv,ind);
    const auto vy_y = smart_deriv<0,1,0>(a, v.yv,ind);
    const auto vz_y = smart_deriv<0,0,1>(a, v.yv,ind);

    const auto vx_z = smart_deriv<1,0,0>(a, v.zv,ind);
    const auto vy_z = smart_deriv<0,1,0>(a, v.zv,ind);
    const auto vz_z = smart_deriv<0,0,1>(a, v.zv,ind);


    const auto vxx_x = smart_deriv<2,0,0>(a, v.xv,ind);
    const auto vyy_y = smart_deriv<0,2,0>(a, v.yv,ind);
    const auto vzz_z = smart_deriv<0,0,2>(a, v.zv,ind);


    const auto vxy_x = smart_deriv<1,1,0>(a, v.xv,ind);

    const auto vxz_x = smart_deriv<1,0,1>(a, v.xv,ind);

    const auto vxy_y = smart_deriv<1,1,0>(a, v.yv,ind);
    const auto vyz_y = smart_deriv<0,1,1>(a, v.yv,ind);

    const auto vxz_z = smart_deriv<1,0,1>(a, v.zv,ind);
    const auto vyz_z = smart_deriv<0,1,1>(a, v.zv,ind);

    const auto term1 = vx_x*vx_x + vy_y*vy_y + vz_z*vz_z;
    const auto term2 = 2*vx_y*vy_x + 2*vx_z*vz_x + 2*vy_z*vz_y;
    const auto term3 = s_v.x() * ( vxx_x + vxy_y + vxz_z );
    const auto term4 = s_v.y() * ( vxy_x + vyy_y + vyz_z );
    const auto term5 = s_v.z() * ( vxz_x + vyz_y + vzz_z );


    return term1 + term2 + term3 + term4 + term5;
}





//D(L(v_{ind}))
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline double divergence_laplacian_old(const big_vec_v& v, const unsigned ind) noexcept {
    const auto vxxx = smart_deriv_old<3,0,0>(v.xv,ind);
    const auto vxyy = smart_deriv_old<1,2,0>(v.xv,ind);
    const auto vxzz = smart_deriv_old<1,0,2>(v.xv,ind);

    const auto vyxx = smart_deriv_old<2,1,0>(v.yv,ind);
    const auto vyyy = smart_deriv_old<0,3,0>(v.yv,ind);
    const auto vyzz = smart_deriv_old<0,1,2>(v.yv,ind);

    const auto vzxx = smart_deriv_old<2,0,1>(v.zv,ind);
    const auto vzyy = smart_deriv_old<0,2,1>(v.zv,ind);
    const auto vzzz = smart_deriv_old<0,0,3>(v.zv,ind);

    return vxxx+ vxyy + vxzz + vyxx + vyyy + vyzz + vzxx + vzyy + vzzz;
}


#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline double divergence_laplacian(const axes&a, const big_vec_v& v, const unsigned ind) noexcept {
    /*const auto vxxx = smart_deriv<3,0,0>(a, v.xv,ind);
    const auto vxyy = smart_deriv<1,2,0>(a, v.xv,ind);
    const auto vxzz = smart_deriv<1,0,2>(a, v.xv,ind);

    const auto vyxx = smart_deriv<2,1,0>(a, v.yv,ind);
    const auto vyyy = smart_deriv<0,3,0>(a, v.yv,ind);
    const auto vyzz = smart_deriv<0,1,2>(a, v.yv,ind);

    const auto vzxx = smart_deriv<2,0,1>(a, v.zv,ind);
    const auto vzyy = smart_deriv<0,2,1>(a, v.zv,ind);
    const auto vzzz = smart_deriv<0,0,3>(a, v.zv,ind);

    return vxxx+ vxyy + vxzz + vyxx + vyyy + vyzz + vzxx + vzyy + vzzz;*/
    return divergence_laplacian_old(v, ind);
}




#endif //CODE_CALC_HPP
