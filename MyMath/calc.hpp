//
// Created by jacob on 11/8/21.
//

#ifndef CODE_CALC_HPP
#define CODE_CALC_HPP

#include "vec3.hpp"
#include "big_vec.hpp"  //to removed later




//nx,ny,nz denote the derivative. e.g. <1,1,0> corresponds to d^2/dxdy
//i,j,k corresponds to where the derivative is to take place
template<unsigned nx, unsigned ny, unsigned nz, typename T>
auto smart_deriv(const T& v, const unsigned i, const unsigned j, const unsigned k) {
    static_assert( (nx < 4) && (ny < 4) && (nz < 4), "Forth derivatives and higher are not supported");
    static_assert( (nx+ny+nz) < 4, "Mixed derivatives higher than third order are not supported" );
    static_assert( !( (nx==1) && (ny==1) && (nz==1) ), "d^3/dxdydz is not currently supported" );
    static_assert( (nx+ny+nz) > 0, "Need at least 1 dimension for a derivative" );

    constexpr auto m = std::max( {nx, ny, nz} );
    constexpr auto s = nx + ny + nz;

    if constexpr (s==1) {   //first order derivative
        if constexpr(nx==1) {   //d/dx
            if (!v.has_left(i,j,k)) {    //at left boundary --- forward difference needed
                return forward_difference_1st<0>(v, i, j, k, v.dx, v.dy, v.dz);
            } else if (!v.has_right(i,j,k))  {   //at right boundary --- need backward difference
                return backward_difference_1st<0>(v, i, j, k, v.dx, v.dy, v.dz);
            } else {    //nowhere special --- can just use central difference
                return central_difference_1st<0>(v, i, j, k, v.dx, v.dy, v.dz);
            }
        } else if constexpr(ny==1) {    //d/dy
            if (!v.has_down(i,j,k)) {    //at bottom boundary --- forward difference needed
                return forward_difference_1st<1>(v, i, j, k, v.dx, v.dy, v.dz);
            } else if (!v.has_up(i,j,k))  {   //at top boundary --- need backward difference
                return backward_difference_1st<1>(v, i, j, k, v.dx, v.dy, v.dz);
            } else {    //nowhere special --- can just use central difference
                return central_difference_1st<1>(v, i, j, k, v.dx, v.dy, v.dz);
            }
        } else if constexpr(nz==1) {    //d/dz
            if (!v.has_front(i,j,k)) {    //at front boundary --- forward difference needed
                return forward_difference_1st<2>(v, i, j, k, v.dx, v.dy, v.dz);
            } else if (!v.has_back(i,j,k))  {   //at back boundary --- need backward difference
                return backward_difference_1st<2>(v, i, j, k, v.dx, v.dy, v.dz);
            } else {    //nowhere special --- can just use central difference
                return central_difference_1st<2>(v, i, j, k, v.dx, v.dy, v.dz);
            }
        }
    }
    else if constexpr(s==2) { //second order derivative
        if constexpr(m==2) {    //pure derivative
            if constexpr(nx==2) {   //d^2/dx^2
                if (!v.has_left(i,j,k)) {    //at left boundary --- forward difference needed
                    return forward_difference_2nd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                } else if (!v.has_right(i,j,k))  {   //at right boundary --- need backward difference
                    return backward_difference_2nd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                } else {    //nowhere special --- can just use central difference
                    //return forward_difference_2nd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                    return central_difference_2nd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                }
            } else if constexpr(ny==2) {    //d^2/dy^2
                if (!v.has_down(i,j,k)) {    //at bottom boundary --- forward difference needed
                    return forward_difference_2nd<1>(v, i, j, k, v.dx, v.dy, v.dz);
                } else if (!v.has_up(i,j,k))  {   //at top boundary --- need backward difference
                    return backward_difference_2nd<1>(v, i, j, k, v.dx, v.dy, v.dz);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_2nd<1>(v, i, j, k, v.dx, v.dy, v.dz);
                }
            } else if constexpr(nz==2) {    //d^2/dz^2
                if (!v.has_front(i,j,k)) {    //at front boundary --- forward difference needed
                    return forward_difference_2nd<2>(v, i, j, k, v.dx, v.dy, v.dz);
                } else if (!v.has_back(i,j,k))  {   //at back boundary --- need backward difference
                    return backward_difference_2nd<2>(v, i, j, k, v.dx, v.dy, v.dz);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_2nd<2>(v, i, j, k, v.dx, v.dy, v.dz);
                }
            }
        }
        else {    //mixed derivative
            if constexpr( (nx==1) && (ny==1)) { //d^2/dxdy
                if (!v.has_left(i,j,k)) {
                    if (!v.has_down(i,j,k)) { //forward forward
                        return forward_difference_2nd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //forward backward
                        return forward_backward_difference_2nd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return central_forward_difference_2nd_mixed<1,0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_right(i,j,k)) {
                    if (!v.has_down(i,j,k)) { //backward forward
                        return forward_backward_difference_2nd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //backward backward
                        return backward_difference_2nd_mixed<0,1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return central_backward_difference_2nd_mixed<1,0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_down(i,j,k)) { //central forward
                        return central_forward_difference_2nd_mixed<0,1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //central backward
                        return central_backward_difference_2nd_mixed<0,1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_2nd_mixed<0,1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            } else if constexpr( (nx==1) && (nz==1)) {  //d^2/dxdz
                if (!v.has_left(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //forward forward
                        return forward_difference_2nd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //forward backward
                        return forward_backward_difference_2nd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return central_forward_difference_2nd_mixed<2,0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_right(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //backward forward
                        return forward_backward_difference_2nd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //backward backward
                        return backward_difference_2nd_mixed<0,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return central_backward_difference_2nd_mixed<2,0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_front(i,j,k)) { //central forward
                        return central_forward_difference_2nd_mixed<0,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //central backward
                        return central_backward_difference_2nd_mixed<0,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_2nd_mixed<0,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            } else if constexpr( (ny==1) && (nz==1)) {  //d^2/dydz
                if (!v.has_down(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //forward forward
                        return forward_difference_2nd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //forward backward
                        return forward_backward_difference_2nd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return central_forward_difference_2nd_mixed<2,1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_up(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //backward forward
                        return forward_backward_difference_2nd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //backward backward
                        return backward_difference_2nd_mixed<1,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return central_backward_difference_2nd_mixed<2,1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_front(i,j,k)) { //central forward
                        return central_forward_difference_2nd_mixed<1,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //central backward
                        return central_backward_difference_2nd_mixed<1,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_2nd_mixed<1,2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }
        }

    }
    else if constexpr(s==3) { //third order derivative
        if constexpr(m==3) {    //pure derivative
            if constexpr(nx==3) {   //d^3/dx^3
                if (!v.has_2left(i,j,k)) {    //at left boundary --- forward difference needed
                    return forward_difference_3rd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                } else if (!v.has_2right(i,j,k))  {   //at right boundary --- need backward difference
                    return backward_difference_3rd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_3rd<0>(v, i, j, k, v.dx, v.dy, v.dz);
                }
            } else if constexpr(ny==3) {    //d^3/dy^3
                if (!v.has_2down(i,j,k)) {    //at bottom boundary --- forward difference needed
                    return forward_difference_3rd<1>(v, i, j, k, v.dx, v.dy, v.dz);
                } else if (!v.has_2up(i,j,k))  {   //at top boundary --- need backward difference
                    return backward_difference_3rd<1>(v, i, j, k, v.dx, v.dy, v.dz);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_3rd<1>(v, i, j, k, v.dx, v.dy, v.dz);
                }
            } else if constexpr(nz==3) {    //d^3/dz^3
                if (!v.has_2front(i,j,k)) {    //at front boundary --- forward difference needed
                    return forward_difference_3rd<2>(v, i, j, k, v.dx, v.dy, v.dz);
                } else if (!v.has_2back(i,j,k))  {   //at back boundary --- need backward difference
                    return backward_difference_3rd<2>(v, i, j, k, v.dx, v.dy, v.dz);
                } else {    //nowhere special --- can just use central difference
                    return central_difference_3rd<2>(v, i, j, k, v.dx, v.dy, v.dz);
                }
            }
        } else {    //mixed derivative
            if constexpr( (nx==2) && (ny==1) ) {    //d^3/dx^2dy
                if (!v.has_left(i,j,k)) {
                    if (!v.has_down(i,j,k)) { //forward forward
                        return forward_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //forward backward
                        return forward_backward_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return forward_central_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_right(i,j,k)) {
                    if (!v.has_down(i,j,k)) { //backward forward
                        return backward_forward_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //backward backward
                        return backward_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return backward_central_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_down(i,j,k)) { //central forward
                        return central_forward_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //central backward
                        return central_backward_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_3rd_mixed<0, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }
            else if constexpr( (nx==1) && (ny==2) ) { //d^3/dy^2dx
                if (!v.has_left(i,j,k)) {
                    if (!v.has_down(i,j,k)) { //forward forward
                        return forward_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //forward backward
                        return backward_forward_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return central_forward_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_right(i,j,k)) {
                    if (!v.has_down(i,j,k)) { //backward forward
                        return forward_backward_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //backward backward
                        return backward_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return central_backward_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_down(i,j,k)) { //central forward
                        return forward_central_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_up(i,j,k)) {    //central backward
                        return backward_central_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_3rd_mixed<1, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }
            else if constexpr( (nx==2) && (nz==1) ) { //d^3/dx^2dz
                if (!v.has_left(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //forward forward
                        return forward_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //forward backward
                        return forward_backward_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return forward_central_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_right(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //backward forward
                        return backward_forward_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //backward backward
                        return backward_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return backward_central_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_front(i,j,k)) { //central forward
                        return central_forward_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //central backward
                        return central_backward_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_3rd_mixed<0, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }
            else if constexpr( (nx==1) && (nz==2) ) { //d^3/dz^2dx
                if (!v.has_left(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //forward forward
                        return forward_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //forward backward
                        return backward_forward_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return central_forward_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_right(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //backward forward
                        return forward_backward_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //backward backward
                        return backward_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return central_backward_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_front(i,j,k)) { //central forward
                        return forward_central_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //central backward
                        return backward_central_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_3rd_mixed<2, 0>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }
            else if constexpr( (nz==2) && (ny==1) ) { //d^3/dz^2dy
                if (!v.has_down(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //forward forward
                        return forward_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //forward backward
                        return backward_forward_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return central_forward_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_up(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //backward forward
                        return forward_backward_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //backward backward
                        return backward_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return central_backward_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_front(i,j,k)) { //central forward
                        return forward_central_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //central backward
                        return backward_central_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_3rd_mixed<2, 1>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }
            else if constexpr( (nz==1) && (ny==2) ) { //d^3/dy^2dz
                if (!v.has_down(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //forward forward
                        return forward_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //forward backward
                        return forward_backward_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //forward central
                        return forward_central_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else if (!v.has_up(i,j,k)) {
                    if (!v.has_front(i,j,k)) { //backward forward
                        return backward_forward_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //backward backward
                        return backward_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //backward central
                        return backward_central_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                } else {
                    if (!v.has_front(i,j,k)) { //central forward
                        return central_forward_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else if (!v.has_back(i,j,k)) {    //central backward
                        return central_backward_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    } else {    //central central
                        return central_difference_3rd_mixed<1, 2>(v, i, j, k, v.dx, v.dy, v.dz);
                    }
                }
            }

        }
    }
}


//H(v_{i,j,k})
template <unsigned N, unsigned M, unsigned P>
vec3 advection(const big_vec<N,M,P, vec3>& v, const unsigned i, const unsigned j, const unsigned k) {
    const auto s_v = v(i,j,k);  //small v
    return s_v.x() * smart_deriv<1,0,0>(v, i,j,k) + s_v.y() * smart_deriv<0,1,0>(v, i,j,k) + s_v.z() * smart_deriv<0,0,1>(v, i,j,k) ;
}

//L(v_{i,j,k})
template <typename T>
auto laplacian(const T& v, const unsigned i, const unsigned j, const unsigned k) {
    return smart_deriv<2,0,0>(v, i,j,k) + smart_deriv<0,2,0>(v, i,j,k) +smart_deriv<0,0,2>(v, i,j,k);
}

//TODO : Test
template <unsigned N, unsigned M, unsigned P>
vec3 gradient(const big_vec<N,M,P, double>& v, const unsigned i, const unsigned j, const unsigned k) {
    return vec3(smart_deriv<1,0,0>(v, i, j, k), smart_deriv<0,1,0>(v, i, j, k), smart_deriv<0,0,1>(v, i, j, k));
}

//D(v_{i,j,k})
template <unsigned N, unsigned M, unsigned P>
double divergence(const big_vec<N,M,P, vec3>& v, const unsigned i, const unsigned j, const unsigned k) {
    return smart_deriv<1,0,0>(v.xv, i, j, k) + smart_deriv<0,1,0>(v.yv, i, j, k) + smart_deriv<0,0,1>(v.zv, i, j, k);
}

//D(H(v_{i,j,k}))
template <unsigned N, unsigned M, unsigned P>
double divergence_advection(const big_vec<N,M,P, vec3>& v, const unsigned i, const unsigned j, const unsigned k) {
    const auto s_v = v(i,j,k);  //small vec

    const auto vx_x = smart_deriv<1,0,0>(v.xv,i,j,k);
    const auto vy_x = smart_deriv<0,1,0>(v.xv,i,j,k);
    const auto vz_x = smart_deriv<0,0,1>(v.xv,i,j,k);

    const auto vx_y = smart_deriv<1,0,0>(v.yv,i,j,k);
    const auto vy_y = smart_deriv<0,1,0>(v.yv,i,j,k);
    const auto vz_y = smart_deriv<0,0,1>(v.yv,i,j,k);

    const auto vx_z = smart_deriv<1,0,0>(v.zv,i,j,k);
    const auto vy_z = smart_deriv<0,1,0>(v.zv,i,j,k);
    const auto vz_z = smart_deriv<0,0,1>(v.zv,i,j,k);



    const auto vxx_x = smart_deriv<2,0,0>(v.xv,i,j,k);

    const auto vyy_y = smart_deriv<0,2,0>(v.yv,i,j,k);

    const auto vzz_z = smart_deriv<0,0,2>(v.zv,i,j,k);




    const auto vxy_x = smart_deriv<1,1,0>(v.xv,i,j,k);
    const auto vxz_x = smart_deriv<1,0,1>(v.xv,i,j,k);

    const auto vxy_y = smart_deriv<1,1,0>(v.yv,i,j,k);
    const auto vyz_y = smart_deriv<0,1,1>(v.yv,i,j,k);

    const auto vxz_z = smart_deriv<1,0,1>(v.zv,i,j,k);
    const auto vyz_z = smart_deriv<0,1,1>(v.zv,i,j,k);



    const auto term1 = vx_x*vx_x + vy_y*vy_y + vz_z*vz_z;
    const auto term2 = 2*vx_y*vy_x + 2*vx_z*vz_x + 2*vy_z*vz_y;
    const auto term3 = s_v.x() * ( vxx_x + vxy_y + vxz_z );
    const auto term4 = s_v.y() * ( vxy_x + vyy_y + vyz_z );
    const auto term5 = s_v.z() * ( vxz_x + vyz_y + vzz_z );


    return term1 + term2 + term3 + term4 + term5;
}

//D(L(v_{i,j,k}))
template <unsigned N, unsigned M, unsigned P>
double divergence_laplacian(const big_vec<N,M,P, vec3>& v, const unsigned i, const unsigned j, const unsigned k) {
    const auto vxxx = smart_deriv<3,0,0>(v.xv,i,j,k);
    const auto vxyy = smart_deriv<1,2,0>(v.xv,i,j,k);
    const auto vxzz = smart_deriv<1,0,2>(v.xv,i,j,k);

    const auto vyxx = smart_deriv<2,1,0>(v.yv,i,j,k);
    const auto vyyy = smart_deriv<0,3,0>(v.yv,i,j,k);
    const auto vyzz = smart_deriv<0,1,2>(v.yv,i,j,k);

    const auto vzxx = smart_deriv<2,0,1>(v.zv,i,j,k);
    const auto vzyy = smart_deriv<0,2,1>(v.zv,i,j,k);
    const auto vzzz = smart_deriv<0,0,3>(v.zv,i,j,k);

    return vxxx+ vxyy + vxzz + vyxx + vyyy + vyzz + vzxx + vzyy + vzzz;
}




#endif //CODE_CALC_HPP
