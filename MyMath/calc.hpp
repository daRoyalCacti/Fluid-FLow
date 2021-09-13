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
    static_assert( !( (nx==1) && (ny==1) && (nz==1) ), "d^3/dxdydz is not currently supported" );
    static_assert( (nx+ny+nz) > 0, "Need at least 1 dimension for a derivative" );

    constexpr auto m = std::max( {nx, ny, nz} );
    constexpr auto s = nx + ny + nz;

    if constexpr (s==1) {   //first order derivative
        if constexpr(nx==1) {   //d/dx
            if (!v.has_left(ind)) {    //at left boundary --- forward difference needed
                return forward_difference_1st<0>(v, ind);
            } else if (!v.has_right(ind))  {   //at right boundary --- need backward difference
                return backward_difference_1st<0>(v, ind);
            } else {    //nowhere special --- can just use central difference
                return central_difference_1st<0>(v, ind);
            }
        } else if constexpr(ny==1) {    //d/dy
            if (!v.has_down(ind)) {    //at bottom boundary --- forward difference needed
                return forward_difference_1st<1>(v, ind);
            } else if (!v.has_up(ind))  {   //at top boundary --- need backward difference
                return backward_difference_1st<1>(v, ind);
            } else {    //nowhere special --- can just use central difference
                return central_difference_1st<1>(v, ind);
            }
        } else if constexpr(nz==1) {    //d/dz
            if (!v.has_front(ind)) {    //at front boundary --- forward difference needed
                return forward_difference_1st<2>(v, ind);
            } else if (!v.has_back(ind))  {   //at back boundary --- need backward difference
                return backward_difference_1st<2>(v, ind);
            } else {    //nowhere special --- can just use central difference
                return central_difference_1st<2>(v, ind);
            }
        }
    }
    else if constexpr(s==2) { //second order derivative
        if constexpr(m==2) {    //pure derivative
            if constexpr(nx==2) {   //d^2/dx^2
                if (!v.has_left(ind)) {    //at left boundary --- forward difference needed
                    return forward_difference_2nd<0>(v, ind);
                } else if (!v.has_right(ind))  {   //at right boundary --- need backward difference
                    return backward_difference_2nd<0>(v, ind);
                } else {    //nowhere special --- can just use central difference
                    //return forward_difference_2nd<0>(v, ind);
                    return central_difference_2nd<0>(v, ind);
                }
            } else if constexpr(ny==2) {    //d^2/dy^2
                if (!v.has_down(ind)) {    //at bottom boundary --- forward difference needed
                    return forward_difference_2nd<1>(v, ind);
                } else if (!v.has_up(ind))  {   //at top boundary --- need backward difference
                    return backward_difference_2nd<1>(v, ind);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_2nd<1>(v, ind);
                }
            } else if constexpr(nz==2) {    //d^2/dz^2
                if (!v.has_front(ind)) {    //at front boundary --- forward difference needed
                    return forward_difference_2nd<2>(v, ind);
                } else if (!v.has_back(ind))  {   //at back boundary --- need backward difference
                    return backward_difference_2nd<2>(v, ind);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_2nd<2>(v, ind);
                }
            }
        }
        else {    //mixed derivative
            if constexpr( (nx==1) && (ny==1)) { //d^2/dxdy
                if (!v.has_left(ind)) {
                    if (!v.has_down(ind)) { //forward forward
                        return forward_difference_2nd_mixed<0, 1>(v, ind);
                    } else if (!v.has_up(ind)) {    //forward backward
                        return forward_backward_difference_2nd_mixed<0, 1>(v, ind);
                    } else {    //forward central
                        return central_forward_difference_2nd_mixed<1,0>(v, ind);
                    }
                } else if (!v.has_right(ind)) {
                    if (!v.has_down(ind)) { //backward forward
                        return forward_backward_difference_2nd_mixed<1, 0>(v, ind);
                    } else if (!v.has_up(ind)) {    //backward backward
                        return backward_difference_2nd_mixed<0,1>(v, ind);
                    } else {    //backward central
                        return central_backward_difference_2nd_mixed<1,0>(v, ind);
                    }
                } else {
                    if (!v.has_down(ind)) { //central forward
                        return central_forward_difference_2nd_mixed<0,1>(v, ind);
                    } else if (!v.has_up(ind)) {    //central backward
                        return central_backward_difference_2nd_mixed<0,1>(v, ind);
                    } else {    //central central
                        return central_difference_2nd_mixed<0,1>(v, ind);
                    }
                }
            } else if constexpr( (nx==1) && (nz==1)) {  //d^2/dxdz
                if (!v.has_left(ind)) {
                    if (!v.has_front(ind)) { //forward forward
                        return forward_difference_2nd_mixed<0, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //forward backward
                        return forward_backward_difference_2nd_mixed<0, 2>(v, ind);
                    } else {    //forward central
                        return central_forward_difference_2nd_mixed<2,0>(v, ind);
                    }
                } else if (!v.has_right(ind)) {
                    if (!v.has_front(ind)) { //backward forward
                        return forward_backward_difference_2nd_mixed<2, 0>(v, ind);
                    } else if (!v.has_back(ind)) {    //backward backward
                        return backward_difference_2nd_mixed<0,2>(v, ind);
                    } else {    //backward central
                        return central_backward_difference_2nd_mixed<2,0>(v, ind);
                    }
                } else {
                    if (!v.has_front(ind)) { //central forward
                        return central_forward_difference_2nd_mixed<0,2>(v, ind);
                    } else if (!v.has_back(ind)) {    //central backward
                        return central_backward_difference_2nd_mixed<0,2>(v, ind);
                    } else {    //central central
                        return central_difference_2nd_mixed<0,2>(v, ind);
                    }
                }
            } else if constexpr( (ny==1) && (nz==1)) {  //d^2/dydz
                if (!v.has_down(ind)) {
                    if (!v.has_front(ind)) { //forward forward
                        return forward_difference_2nd_mixed<1, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //forward backward
                        return forward_backward_difference_2nd_mixed<1, 2>(v, ind);
                    } else {    //forward central
                        return central_forward_difference_2nd_mixed<2,1>(v, ind);
                    }
                } else if (!v.has_up(ind)) {
                    if (!v.has_front(ind)) { //backward forward
                        return forward_backward_difference_2nd_mixed<2, 1>(v, ind);
                    } else if (!v.has_back(ind)) {    //backward backward
                        return backward_difference_2nd_mixed<1,2>(v, ind);
                    } else {    //backward central
                        return central_backward_difference_2nd_mixed<2,1>(v, ind);
                    }
                } else {
                    if (!v.has_front(ind)) { //central forward
                        return central_forward_difference_2nd_mixed<1,2>(v, ind);
                    } else if (!v.has_back(ind)) {    //central backward
                        return central_backward_difference_2nd_mixed<1,2>(v, ind);
                    } else {    //central central
                        return central_difference_2nd_mixed<1,2>(v, ind);
                    }
                }
            }
        }

    }
    else if constexpr(s==3) { //third order derivative
        if constexpr(m==3) {    //pure derivative
            if constexpr(nx==3) {   //d^3/dx^3
                if (!v.has_2left(ind)) {    //at left boundary --- forward difference needed
                    return forward_difference_3rd<0>(v, ind);
                } else if (!v.has_2right(ind))  {   //at right boundary --- need backward difference
                    return backward_difference_3rd<0>(v, ind);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_3rd<0>(v, ind);
                }
            } else if constexpr(ny==3) {    //d^3/dy^3
                if (!v.has_2down(ind)) {    //at bottom boundary --- forward difference needed
                    return forward_difference_3rd<1>(v, ind);
                } else if (!v.has_2up(ind))  {   //at top boundary --- need backward difference
                    return backward_difference_3rd<1>(v, ind);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_3rd<1>(v, ind);
                }
            } else if constexpr(nz==3) {    //d^3/dz^3
                if (!v.has_2front(ind)) {    //at front boundary --- forward difference needed
                    return forward_difference_3rd<2>(v, ind);
                } else if (!v.has_2back(ind))  {   //at back boundary --- need backward difference
                    return backward_difference_3rd<2>(v, ind);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_3rd<2>(v, ind);
                }
            }
        } else {    //mixed derivative
            if constexpr( (nx==2) && (ny==1) ) {    //d^3/dx^2dy
                if (!v.has_left(ind)) {
                    if (!v.has_down(ind)) { //forward forward
                        return forward_difference_3rd_mixed<0, 1>(v, ind);
                    } else if (!v.has_up(ind)) {    //forward backward
                        return forward_backward_difference_3rd_mixed<0, 1>(v, ind);
                    } else {    //forward central
                        return forward_central_difference_3rd_mixed<0, 1>(v, ind);
                    }
                } else if (!v.has_right(ind)) {
                    if (!v.has_down(ind)) { //backward forward
                        return backward_forward_difference_3rd_mixed<0, 1>(v, ind);
                    } else if (!v.has_up(ind)) {    //backward backward
                        return backward_difference_3rd_mixed<0, 1>(v, ind);
                    } else {    //backward central
                        return backward_central_difference_3rd_mixed<0, 1>(v, ind);
                    }
                } else {
                    if (!v.has_down(ind)) { //central forward
                        return central_forward_difference_3rd_mixed<0, 1>(v, ind);
                    } else if (!v.has_up(ind)) {    //central backward
                        return central_backward_difference_3rd_mixed<0, 1>(v, ind);
                    } else {    //central central
                        return central_difference_3rd_mixed<0, 1>(v, ind);
                    }
                }
            }
            else if constexpr( (nx==1) && (ny==2) ) { //d^3/dy^2dx
                if (!v.has_left(ind)) {
                    if (!v.has_down(ind)) { //forward forward
                        return forward_difference_3rd_mixed<1, 0>(v, ind);
                    } else if (!v.has_up(ind)) {    //forward backward
                        return backward_forward_difference_3rd_mixed<1, 0>(v, ind);
                    } else {    //forward central
                        return central_forward_difference_3rd_mixed<1, 0>(v, ind);
                    }
                } else if (!v.has_right(ind)) {
                    if (!v.has_down(ind)) { //backward forward
                        return forward_backward_difference_3rd_mixed<1, 0>(v, ind);
                    } else if (!v.has_up(ind)) {    //backward backward
                        return backward_difference_3rd_mixed<1, 0>(v, ind);
                    } else {    //backward central
                        return central_backward_difference_3rd_mixed<1, 0>(v, ind);
                    }
                } else {
                    if (!v.has_down(ind)) { //central forward
                        return forward_central_difference_3rd_mixed<1, 0>(v, ind);
                    } else if (!v.has_up(ind)) {    //central backward
                        return backward_central_difference_3rd_mixed<1, 0>(v, ind);
                    } else {    //central central
                        return central_difference_3rd_mixed<1, 0>(v, ind);
                    }
                }
            }
            else if constexpr( (nx==2) && (nz==1) ) { //d^3/dx^2dz
                if (!v.has_left(ind)) {
                    if (!v.has_front(ind)) { //forward forward
                        return forward_difference_3rd_mixed<0, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //forward backward
                        return forward_backward_difference_3rd_mixed<0, 2>(v, ind);
                    } else {    //forward central
                        return forward_central_difference_3rd_mixed<0, 2>(v, ind);
                    }
                } else if (!v.has_right(ind)) {
                    if (!v.has_front(ind)) { //backward forward
                        return backward_forward_difference_3rd_mixed<0, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //backward backward
                        return backward_difference_3rd_mixed<0, 2>(v, ind);
                    } else {    //backward central
                        return backward_central_difference_3rd_mixed<0, 2>(v, ind);
                    }
                } else {
                    if (!v.has_front(ind)) { //central forward
                        return central_forward_difference_3rd_mixed<0, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //central backward
                        return central_backward_difference_3rd_mixed<0, 2>(v, ind);
                    } else {    //central central
                        return central_difference_3rd_mixed<0, 2>(v, ind);
                    }
                }
            }
            else if constexpr( (nx==1) && (nz==2) ) { //d^3/dz^2dx
                if (!v.has_left(ind)) {
                    if (!v.has_front(ind)) { //forward forward
                        return forward_difference_3rd_mixed<2, 0>(v, ind);
                    } else if (!v.has_back(ind)) {    //forward backward
                        return backward_forward_difference_3rd_mixed<2, 0>(v, ind);
                    } else {    //forward central
                        return central_forward_difference_3rd_mixed<2, 0>(v, ind);
                    }
                } else if (!v.has_right(ind)) {
                    if (!v.has_front(ind)) { //backward forward
                        return forward_backward_difference_3rd_mixed<2, 0>(v, ind);
                    } else if (!v.has_back(ind)) {    //backward backward
                        return backward_difference_3rd_mixed<2, 0>(v, ind);
                    } else {    //backward central
                        return central_backward_difference_3rd_mixed<2, 0>(v, ind);
                    }
                } else {
                    if (!v.has_front(ind)) { //central forward
                        return forward_central_difference_3rd_mixed<2, 0>(v, ind);
                    } else if (!v.has_back(ind)) {    //central backward
                        return backward_central_difference_3rd_mixed<2, 0>(v, ind);
                    } else {    //central central
                        return central_difference_3rd_mixed<2, 0>(v, ind);
                    }
                }
            }
            else if constexpr( (nz==2) && (ny==1) ) { //d^3/dz^2dy
                if (!v.has_down(ind)) {
                    if (!v.has_front(ind)) { //forward forward
                        return forward_difference_3rd_mixed<2, 1>(v, ind);
                    } else if (!v.has_back(ind)) {    //forward backward
                        return backward_forward_difference_3rd_mixed<2, 1>(v, ind);
                    } else {    //forward central
                        return central_forward_difference_3rd_mixed<2, 1>(v, ind);
                    }
                } else if (!v.has_up(ind)) {
                    if (!v.has_front(ind)) { //backward forward
                        return forward_backward_difference_3rd_mixed<2, 1>(v, ind);
                    } else if (!v.has_back(ind)) {    //backward backward
                        return backward_difference_3rd_mixed<2, 1>(v, ind);
                    } else {    //backward central
                        return central_backward_difference_3rd_mixed<2, 1>(v, ind);
                    }
                } else {
                    if (!v.has_front(ind)) { //central forward
                        return forward_central_difference_3rd_mixed<2, 1>(v, ind);
                    } else if (!v.has_back(ind)) {    //central backward
                        return backward_central_difference_3rd_mixed<2, 1>(v, ind);
                    } else {    //central central
                        return central_difference_3rd_mixed<2, 1>(v, ind);
                    }
                }
            }
            else if constexpr( (nz==1) && (ny==2) ) { //d^3/dy^2dz
                if (!v.has_down(ind)) {
                    if (!v.has_front(ind)) { //forward forward
                        return forward_difference_3rd_mixed<1, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //forward backward
                        return forward_backward_difference_3rd_mixed<1, 2>(v, ind);
                    } else {    //forward central
                        return forward_central_difference_3rd_mixed<1, 2>(v, ind);
                    }
                } else if (!v.has_up(ind)) {
                    if (!v.has_front(ind)) { //backward forward
                        return backward_forward_difference_3rd_mixed<1, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //backward backward
                        return backward_difference_3rd_mixed<1, 2>(v, ind);
                    } else {    //backward central
                        return backward_central_difference_3rd_mixed<1, 2>(v, ind);
                    }
                } else {
                    if (!v.has_front(ind)) { //central forward
                        return central_forward_difference_3rd_mixed<1, 2>(v, ind);
                    } else if (!v.has_back(ind)) {    //central backward
                        return central_backward_difference_3rd_mixed<1, 2>(v, ind);
                    } else {    //central central
                        return central_difference_3rd_mixed<1, 2>(v, ind);
                    }
                }
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
template <unsigned N, unsigned M, unsigned P>
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
