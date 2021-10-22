//
// Created by jacob on 8/8/21.
//

#ifndef CODE_FINITE_DIFFERENCE_HPP
#define CODE_FINITE_DIFFERENCE_HPP
/*  Axis labelling
 *  Axis = 0 -> x axis
 *  Axis = 1 -> y axis
 *  Axis = 2 -> z axis
 */


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_difference_1st(const T& data, const size_t index) noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data.move(index,1,0,0) - data.move(index,-1,0,0) ) / (2*data.dx(index));
    } else if constexpr(axis == 1) {
        return ( data.move(index,0,1,0) - data.move(index,0,-1,0) ) / (2*data.dy(index));
    } else if constexpr(axis == 2) {
        return ( data.move(index,0,0,1) - data.move(index,0,0,-1) ) / (2*data.dz(index));
    }
}

template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_difference_1st(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( -data.move(index,2,0,0) + 4*data.move(index,1,0,0) - 3*data(index)) / (2*data.dx(index));
    } else if constexpr(axis == 1) {
        return ( -data.move(index,0,2,0) + 4*data.move(index,0,1,0) - 3*data(index)) / (2*data.dy(index));
    } else if constexpr(axis == 2) {
        return ( -data.move(index,0,0,2) + 4*data.move(index,0,0,1) - 3*data(index)) / (2*data.dz(index));
    }
}


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_difference_1st(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( 3*data(index) - 4*data.move(index,-1,0,0) + data.move(index,-2,0,0) ) / (2*data.dx(index));
    } else if constexpr(axis == 1) {
        return ( 3*data(index) - 4*data.move(index,0,-1,0) + data.move(index,0,-2,0) ) / (2*data.dy(index));
    } else if constexpr(axis == 2) {
        return ( 3*data(index) - 4*data.move(index,0,0,-1) + data.move(index,0,0,-2) ) / (2*data.dz(index));
    }
}


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_difference_2nd(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data.move(index,1,0,0) - 2* data(index) + data.move(index,-1,0,0) ) / (data.dx(index)*data.dx(index));
    } else if constexpr(axis == 1) {
        return  ( data.move(index,0,1,0) - 2* data(index) + data.move(index,0,-1,0) ) / (data.dy(index)*data.dy(index));
    } else if constexpr(axis == 2) {
        return  ( data.move(index,0,0,1) - 2* data(index) + data.move(index,0,0,-1) ) / (data.dz(index)*data.dz(index));
    }
}


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_difference_2nd(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (-data.move(index,3,0,0)+4*data.move(index,2,0,0) - 5*data.move(index,1,0,0)+2*data(index)) / (data.dx(index)*data.dx(index));
    } else if constexpr(axis == 1) {
        return (-data.move(index,0,3,0)+4*data.move(index,0,2,0) - 5*data.move(index,0,1,0)+2*data(index)) / (data.dy(index)*data.dy(index));
    } else if constexpr(axis == 2) {
        return (-data.move(index,0,0,3)+4*data.move(index,0,0,2) - 5*data.move(index,0,0,1)+2*data(index)) / (data.dz(index)*data.dz(index));
    }
}


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_difference_2nd(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( 2*data(index) -5*data.move(index,-1,0,0)+4*data.move(index,-2,0,0) - data.move(index,-3,0,0) ) / (data.dx(index)*data.dx(index));
    } else if constexpr(axis == 1) {
        return ( 2*data(index) -5*data.move(index,0,-1,0)+4*data.move(index,0,-2,0) - data.move(index,0,-3,0) ) / (data.dy(index)*data.dy(index));
    } else if constexpr(axis == 2) {
        return ( 2*data(index) -5*data.move(index,0,0,-1)+4*data.move(index,0,0,-2) - data.move(index,0,0,-3) ) / (data.dz(index)*data.dz(index));
    }
}

template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_difference_2nd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");


    if constexpr(axis1 == axis2) {
        return central_difference_2nd<axis1>(data, index);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data.move(index,1,1,0) - data.move(index,-1,1,0) - data.move(index,1,-1,0) + data.move(index,-1,-1,0))  / (4*data.dx(index)*data.dy(index));
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data.move(index,1,0,1) - data.move(index,-1,0,1) - data.move(index,1,0,-1) + data.move(index,-1,0,-1))  / (4*data.dx(index)*data.dz(index));
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data.move(index,0,1,1) - data.move(index,0,-1,1) - data.move(index,0,1,-1) + data.move(index,0,-1,-1))  / (4*data.dy(index)*data.dz(index));
    }
}

template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_difference_2nd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return forward_difference_2nd<axis1>(data, index);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data.move(index,2,2,0) - 4*data.move(index,1,2,0) + 3*data.move(index,0,2,0) - 4*data.move(index,2,1,0)+16*data.move(index,1,1,0)
                    -12*data.move(index,0,1,0) + 3*data.move(index,2,0,0) - 12*data.move(index,1,0,0) + 9*data(index)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data.move(index,2,0,2) - 4*data.move(index,1,0,2) + 3*data.move(index,0,0,2) - 4*data.move(index,2,0,1)+16*data.move(index,1,0,1)
                 -12*data.move(index,0,0,1) + 3*data.move(index,2,0,0) - 12*data.move(index,1,0,0) + 9*data(index)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data.move(index,0,2,2) - 4*data.move(index,0,1,2) + 3*data.move(index,0,0,2) - 4*data.move(index,0,2,1)+16*data.move(index,0,1,1)
                 -12*data.move(index,0,0,1) + 3*data.move(index,0,2,0) - 12*data.move(index,0,1,0) + 9*data(index)) / (4*data.dy(index)*data.dz(index));
    }
}


template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_difference_2nd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return backward_difference_2nd<axis1>(data, index);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data.move(index,-2,-2,0) - 4*data.move(index,-1,-2,0) + 3*data.move(index,0,-2,0) - 4*data.move(index,-2,-1,0)+16*data.move(index,-1,-1,0)
                 -12*data.move(index,0,-1,0) + 3*data.move(index,-2,0,0) - 12*data.move(index,-1,0,0) + 9*data(index)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data.move(index,-2,0,-2) - 4*data.move(index,-1,0,-2) + 3*data.move(index,0,0,-2) - 4*data.move(index,-2,0,-1)+16*data.move(index,-1,0,-1)
                 -12*data.move(index,0,0,-1) + 3*data.move(index,-2,0,0) - 12*data.move(index,-1,0,0) + 9*data(index)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data.move(index,0,-2,-2) - 4*data.move(index,0,-1,-2) + 3*data.move(index,0,0,-2) - 4*data.move(index,0,-2,-1)+16*data.move(index,0,-1,-1)
                 -12*data.move(index,0,0,-1) + 3*data.move(index,0,-2,0) - 12*data.move(index,0,-1,0) + 9*data(index)) / (4*data.dy(index)*data.dz(index));
    }
}



template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_forward_difference_2nd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return (-data.move(index,2,1,0) + 4*data.move(index,1,1,0) - 3*data.move(index,0,1,0) + data.move(index,2,-1,0) - 4*data.move(index,1,-1,0)+3*data.move(index,0,-1,0)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr(axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return (-data.move(index,1,2,0) + 4*data.move(index,1,1,0) - 3*data.move(index,1,0,0) + data.move(index,-1,2,0) - 4*data.move(index,-1,1,0)+3*data.move(index,-1,0,0)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr( axis1 == 2 && axis2 == 0) {//d^2/dxdz
        return (-data.move(index,2,0,1) + 4*data.move(index,1,0,1) - 3*data.move(index,0,0,1) + data.move(index,2,0,-1) - 4*data.move(index,1,0,-1)+3*data.move(index,0,0,-1)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr(axis1 == 0 && axis2 == 2 ) {  //d^2/dxdz
        return (-data.move(index,1,0,2) + 4*data.move(index,1,0,1) - 3*data.move(index,1,0,0) + data.move(index,-1,0,2) - 4*data.move(index,-1,0,1)+3*data.move(index,-1,0,0)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr( axis1 == 2 && axis2 == 1) {//d^2/dydz
        return (-data.move(index,0,2,1) + 4*data.move(index,0,1,1) - 3*data.move(index,0,0,1) + data.move(index,0,2,-1) - 4*data.move(index,0,1,-1)+3*data.move(index,0,0,-1)) / (4*data.dy(index)*data.dz(index));
    } else if constexpr(axis1 == 1 && axis2 == 2) { //d^2/dydz
        return (-data.move(index,0,1,2) + 4*data.move(index,0,1,1) - 3*data.move(index,0,1,0) + data.move(index,0,-1,2) - 4*data.move(index,0,-1,1)+3*data.move(index,0,-1,0)) / (4*data.dy(index)*data.dz(index));
    }
}


template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_backward_difference_2nd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return (data.move(index,-2,1,0) - 4*data.move(index,-1,1,0) + 3*data.move(index,0,1,0) - data.move(index,-2,-1,0) + 4*data.move(index,-1,-1,0)-3*data.move(index,0,-1,0)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr(axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return (data.move(index,1,-2,0) - 4*data.move(index,1,-1,0) + 3*data.move(index,1,0,0) - data.move(index,-1,-2,0) + 4*data.move(index,-1,-1,0)-3*data.move(index,-1,0,0)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr( axis1 == 2 && axis2 == 0) {//d^2/dxdz
        return (data.move(index,-2,0,1) - 4*data.move(index,-1,0,1) + 3*data.move(index,0,0,1) - data.move(index,-2,0,-1) + 4*data.move(index,-1,0,-1)-3*data.move(index,0,0,-1)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr(axis1 == 0 && axis2 == 2 ) {  //d^2/dxdz
        return (data.move(index,1,0,-2) - 4*data.move(index,1,0,-1) + 3*data.move(index,1,0,0) - data.move(index,-1,0,-2) + 4*data.move(index,-1,0,-1)-3*data.move(index,-1,0,0)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr( axis1 == 2 && axis2 == 1) {//d^2/dydz
        return (data.move(index,0,-2,1) - 4*data.move(index,0,-1,1) + 3*data.move(index,0,0,1) - data.move(index,0,-2,-1) + 4*data.move(index,0,-1,-1)-3*data.move(index,0,0,-1)) / (4*data.dy(index)*data.dz(index));
    } else if constexpr(axis1 == 1 && axis2 == 2) { //d^2/dydz
        return (data.move(index,0,1,-2) - 4*data.move(index,0,1,-1) + 3*data.move(index,0,1,0) - data.move(index,0,-1,-2) + 4*data.move(index,0,-1,-1)-3*data.move(index,0,-1,0)) / (4*data.dy(index)*data.dz(index));
    }
}


template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_backward_difference_2nd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return ( -data.move(index,2,-2,0) + 4*data.move(index,1,-2,0) - 3*data.move(index,0,-2,0) + 4*data.move(index,2,-1,0)-16*data.move(index,1,-1,0)
                 +12*data.move(index,0,-1,0) - 3*data.move(index,2,0,0) + 12*data.move(index,1,0,0) - 9*data(index)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return ( -data.move(index,-2,2,0) + 4*data.move(index,-1,2,0) - 3*data.move(index,0,2,0) + 4*data.move(index,-2,1,0)-16*data.move(index,-1,1,0)
                 +12*data.move(index,0,1,0) - 3*data.move(index,-2,0,0) + 12*data.move(index,-1,0,0) - 9*data(index)) / (4*data.dx(index)*data.dy(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^2/dxdz
        return ( -data.move(index,2,0,-2) + 4*data.move(index,1,0,-2) - 3*data.move(index,0,0,-2) + 4*data.move(index,2,0,-1)-16*data.move(index,1,0,-1)
                 +12*data.move(index,0,0,-1) - 3*data.move(index,2,0,0) + 12*data.move(index,1,0,0) - 9*data(index)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^2/dxdz
        return ( -data.move(index,-2,0,2) + 4*data.move(index,-1,0,2) - 3*data.move(index,0,0,2) + 4*data.move(index,-2,0,1)-16*data.move(index,-1,0,1)
                 +12*data.move(index,0,0,1) - 3*data.move(index,-2,0,0) + 12*data.move(index,-1,0,0) - 9*data(index)) / (4*data.dx(index)*data.dz(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^2/dydz
        return ( -data.move(index,0,2,-2) + 4*data.move(index,0,1,-2) - 3*data.move(index,0,0,-2) + 4*data.move(index,0,2,-1)-16*data.move(index,0,1,-1)
                 +12*data.move(index,0,0,-1) - 3*data.move(index,0,2,0) + 12*data.move(index,0,1,0) - 9*data(index)) / (4*data.dy(index)*data.dz(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^2/dydz
        return ( -data.move(index,0,-2,2) + 4*data.move(index,0,-1,2) - 3*data.move(index,0,0,2) + 4*data.move(index,0,-2,1)-16*data.move(index,0,-1,1)
                 +12*data.move(index,0,0,1) - 3*data.move(index,0,-2,0) + 12*data.move(index,0,-1,0) - 9*data(index)) / (4*data.dy(index)*data.dz(index));
    }
}



template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_difference_3rd(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data.move(index,2,0,0) - 2*data.move(index,1,0,0) + 2*data.move(index,-1,0,0) - data.move(index,-2,0,0)) / (2*data.dx(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis == 1) {
        return ( data.move(index,0,2,0) - 2*data.move(index,0,1,0) + 2*data.move(index,0,-1,0) - data.move(index,0,-2,0)) / (2*data.dy(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis == 2) {
        return ( data.move(index,0,0,2) - 2*data.move(index,0,0,1) + 2*data.move(index,0,0,-1) - data.move(index,0,0,-2)) / (2*data.dz(index)*data.dz(index)*data.dz(index));
    }
}


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_difference_3rd(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (-3*data.move(index,4,0,0)+14*data.move(index,3,0,0) - 24*data.move(index,2,0,0)+ 18*data.move(index,1,0,0) - 5*data(index)) / (2*data.dx(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis == 1) {
        return (-3*data.move(index,0,4,0)+14*data.move(index,0,3,0) - 24*data.move(index,0,2,0)+ 18*data.move(index,0,1,0) - 5*data(index)) / (2*data.dy(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis == 2) {
        return (-3*data.move(index,0,0,4)+14*data.move(index,0,0,3) - 24*data.move(index,0,0,2)+ 18*data.move(index,0,0,1) - 5*data(index)) / (2*data.dz(index)*data.dz(index)*data.dz(index));
    }
}


template <unsigned axis, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_difference_3rd(const T& data, const size_t index)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (5*data(index) - 18*data.move(index,-1,0,0) + 24*data.move(index,-2,0,0) - 14*data.move(index,-3,0,0) + 3*data.move(index,-4,0,0)) / (2*data.dx(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis == 1) {
        return (5*data(index) - 18*data.move(index,0,-1,0) + 24*data.move(index,0,-2,0) - 14*data.move(index,0,-3,0) + 3*data.move(index,0,-4,0)) / (2*data.dy(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis == 2) {
        return (5*data(index) - 18*data.move(index,0,0,-1) + 24*data.move(index,0,0,-2) - 14*data.move(index,0,0,-3) + 3*data.move(index,0,0,-4)) / (2*data.dz(index)*data.dz(index)*data.dz(index));
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return central_difference_3rd<axis1>(data, index);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(index,1,1,0) - 2*data.move(index,0,1,0) + data.move(index,-1,1,0) - data.move(index,1,-1,0)+2*data.move(index,0,-1,0)-data.move(index,-1,-1,0)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(index,1,1,0) - 2*data.move(index,1,0,0) + data.move(index,1,-1,0) - data.move(index,-1,1,0)+2*data.move(index,-1,0,0)-data.move(index,-1,-1,0)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(index,1,0,1) - 2*data.move(index,0,0,1) + data.move(index,-1,0,1) - data.move(index,1,0,-1)+2*data.move(index,0,0,-1)-data.move(index,-1,0,-1)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(index,1,0,1) - 2*data.move(index,1,0,0) + data.move(index,1,0,-1) - data.move(index,-1,0,1)+2*data.move(index,-1,0,0)-data.move(index,-1,0,-1)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(index,0,1,1) - 2*data.move(index,0,0,1) + data.move(index,0,-1,1) - data.move(index,0,1,-1)+2*data.move(index,0,0,-1)-data.move(index,0,-1,-1)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(index,0,1,1) - 2*data.move(index,0,1,0) + data.move(index,0,1,-1) - data.move(index,0,-1,1)+2*data.move(index,0,-1,0)-data.move(index,0,-1,-1)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return forward_difference_3rd<axis1>(data, index);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(index,3,2,0) - 4*data.move(index,2,2,0) + 5*data.move(index,1,2,0) - 2*data.move(index,0,2,0) - 4*data.move(index,3,1,0)
                 + 16*data.move(index,2,1,0) - 20*data.move(index,1,1,0) + 8*data.move(index,0,1,0) + 3*data.move(index,3,0,0) - 12*data.move(index,2,0,0)
                 + 15*data.move(index,1,0,0) - 6*data(index)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(index,2,3,0) - 4*data.move(index,2,2,0) + 5*data.move(index,2,1,0) - 2*data.move(index,2,0,0) - 4*data.move(index,1,3,0)
                + 16*data.move(index,1,2,0) - 20*data.move(index,1,1,0) + 8*data.move(index,1,0,0) + 3*data.move(index,0,3,0) - 12*data.move(index,0,2,0)
                + 15*data.move(index,0,1,0) - 6*data(index)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(index,3,0,2) - 4*data.move(index,2,0,2) + 5*data.move(index,1,0,2) - 2*data.move(index,0,0,2) - 4*data.move(index,3,0,1)
                + 16*data.move(index,2,0,1) - 20*data.move(index,1,0,1) + 8*data.move(index,0,0,1) + 3*data.move(index,3,0,0) - 12*data.move(index,2,0,0)
                + 15*data.move(index,1,0,0) - 6*data(index)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(index,2,0,3) - 4*data.move(index,2,0,2) + 5*data.move(index,2,0,1) - 2*data.move(index,2,0,0) - 4*data.move(index,1,0,3)
                + 16*data.move(index,1,0,2) - 20*data.move(index,1,0,1) + 8*data.move(index,1,0,0) + 3*data.move(index,0,0,3) - 12*data.move(index,0,0,2)
                + 15*data.move(index,0,0,1) - 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(index,0,3,2) - 4*data.move(index,0,2,2) + 5*data.move(index,0,1,2) - 2*data.move(index,0,0,2) - 4*data.move(index,0,3,1)
                + 16*data.move(index,0,2,1) - 20*data.move(index,0,1,1) + 8*data.move(index,0,0,1) + 3*data.move(index,0,3,0) - 12*data.move(index,0,2,0)
                + 15*data.move(index,0,1,0) - 6*data(index)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(index,0,2,3) - 4*data.move(index,0,2,2) + 5*data.move(index,0,2,1) - 2*data.move(index,0,2,0) - 4*data.move(index,0,1,3)
                + 16*data.move(index,0,1,2) - 20*data.move(index,0,1,1) + 8*data.move(index,0,1,0) + 3*data.move(index,0,0,3) - 12*data.move(index,0,0,2)
                + 15*data.move(index,0,0,1) - 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return backward_difference_3rd<axis1>(data, index);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(index,-3,-2,0) + 4*data.move(index,-2,-2,0) - 5*data.move(index,-1,-2,0) + 2*data.move(index,0,-2,0) + 4*data.move(index,-3,-1,0)
                - 16*data.move(index,-2,-1,0) + 20*data.move(index,-1,-1,0) - 8*data.move(index,0,-1,0) - 3*data.move(index,-3,0,0) + 12*data.move(index,-2,0,0)
                - 15*data.move(index,-1,0,0) + 6*data(index)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(index,-2,-3,0) + 4*data.move(index,-2,-2,0) - 5*data.move(index,-2,-1,0) + 2*data.move(index,-2,0,0) + 4*data.move(index,-1,-3,0)
                - 16*data.move(index,-1,-2,0) + 20*data.move(index,-1,-1,0) - 8*data.move(index,-1,0,0) - 3*data.move(index,0,-3,0) + 12*data.move(index,0,-2,0)
                - 15*data.move(index,0,-1,0) + 6*data(index)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(index,-3,0,-2) + 4*data.move(index,-2,0,-2) - 5*data.move(index,-1,0,-2) + 2*data.move(index,0,0,-2) + 4*data.move(index,-3,0,-1)
                - 16*data.move(index,-2,0,-1) + 20*data.move(index,-1,0,-1) - 8*data.move(index,0,0,-1) - 3*data.move(index,-3,0,0) + 12*data.move(index,-2,0,0)
                - 15*data.move(index,-1,0,0) + 6*data(index)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(index,-2,0,-3) + 4*data.move(index,-2,0,-2) - 5*data.move(index,-2,0,-1) + 2*data.move(index,-2,0,0) + 4*data.move(index,-1,0,-3)
                - 16*data.move(index,-1,0,-2) + 20*data.move(index,-1,0,-1) - 8*data.move(index,-1,0,0) - 3*data.move(index,0,0,-3) + 12*data.move(index,0,0,-2)
                - 15*data.move(index,0,0,-1) + 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(index,0,-3,-2) + 4*data.move(index,0,-2,-2) - 5*data.move(index,0,-1,-2) + 2*data.move(index,0,0,-2) + 4*data.move(index,0,-3,-1)
                - 16*data.move(index,0,-2,-1) + 20*data.move(index,0,-1,-1) - 8*data.move(index,0,0,-1) - 3*data.move(index,0,-3,0) + 12*data.move(index,0,-2,0)
                - 15*data.move(index,0,-1,0) + 6*data(index)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(index,0,-2,-3) + 4*data.move(index,0,-2,-2) - 5*data.move(index,0,-2,-1) + 2*data.move(index,0,-2,0) + 4*data.move(index,0,-1,-3)
                - 16*data.move(index,0,-1,-2) + 20*data.move(index,0,-1,-1) - 8*data.move(index,0,-1,0) - 3*data.move(index,0,0,-3) + 12*data.move(index,0,0,-2)
                - 15*data.move(index,0,0,-1) + 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_forward_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(index,1,2,0) + 2*data.move(index,0,2,0) - data.move(index,-1,2,0) + 4*data.move(index,1,1,0) - 8*data.move(index,0,1,0)
                 + 4*data.move(index,-1,1,0) - 3*data.move(index,1,0,0) + 6*data(index) - 3*data.move(index,-1,0,0)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(index,2,1,0) + 2*data.move(index,2,0,0) - data.move(index,2,-1,0) + 4*data.move(index,1,1,0) - 8*data.move(index,1,0,0)
                + 4*data.move(index,1,-1,0) - 3*data.move(index,0,1,0) + 6*data(index) - 3*data.move(index,0,-1,0)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(index,1,0,2) + 2*data.move(index,0,0,2) - data.move(index,-1,0,2) + 4*data.move(index,1,0,1) - 8*data.move(index,0,0,1)
                + 4*data.move(index,-1,0,1) - 3*data.move(index,1,0,0) + 6*data(index) - 3*data.move(index,-1,0,0)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(index,2,0,1) + 2*data.move(index,2,0,0) - data.move(index,2,0,-1) + 4*data.move(index,1,0,1) - 8*data.move(index,1,0,0)
                + 4*data.move(index,1,0,-1) - 3*data.move(index,0,0,1) + 6*data(index) - 3*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(index,0,1,2) + 2*data.move(index,0,0,2) - data.move(index,0,-1,2) + 4*data.move(index,0,1,1) - 8*data.move(index,0,0,1)
                + 4*data.move(index,0,-1,1) - 3*data.move(index,0,1,0) + 6*data(index) - 3*data.move(index,0,-1,0)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(index,0,2,1) + 2*data.move(index,0,2,0) - data.move(index,0,2,-1) + 4*data.move(index,0,1,1) - 8*data.move(index,0,1,0)
                + 4*data.move(index,0,1,-1) - 3*data.move(index,0,0,1) + 6*data(index) - 3*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto central_backward_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(index,1,-2,0) - 2*data.move(index,0,-2,0) + data.move(index,-1,-2,0) - 4*data.move(index,1,-1,0) + 8*data.move(index,0,-1,0)
                - 4*data.move(index,-1,-1,0) + 3*data.move(index,1,0,0) - 6*data(index) + 3*data.move(index,-1,0,0)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(index,-2,1,0) - 2*data.move(index,-2,0,0) + data.move(index,-2,-1,0) - 4*data.move(index,-1,1,0) + 8*data.move(index,-1,0,0)
                - 4*data.move(index,-1,-1,0) + 3*data.move(index,0,1,0) - 6*data(index) + 3*data.move(index,0,-1,0)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(index,1,0,-2) - 2*data.move(index,0,0,-2) + data.move(index,-1,0,-2) - 4*data.move(index,1,0,-1) + 8*data.move(index,0,0,-1)
                - 4*data.move(index,-1,0,-1) + 3*data.move(index,1,0,0) - 6*data(index) + 3*data.move(index,-1,0,0)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(index,-2,0,1) - 2*data.move(index,-2,0,0) + data.move(index,-2,0,-1) - 4*data.move(index,-1,0,1) + 8*data.move(index,-1,0,0)
                - 4*data.move(index,-1,0,-1) + 3*data.move(index,0,0,1) - 6*data(index) + 3*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(index,0,1,-2) - 2*data.move(index,0,0,-2) + data.move(index,0,-1,-2) - 4*data.move(index,0,1,-1) + 8*data.move(index,0,0,-1)
                - 4*data.move(index,0,-1,-1) + 3*data.move(index,0,1,0) - 6*data(index) + 3*data.move(index,0,-1,0)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(index,0,-2,1) - 2*data.move(index,0,-2,0) + data.move(index,0,-2,-1) - 4*data.move(index,0,-1,1) + 8*data.move(index,0,-1,0)
                - 4*data.move(index,0,-1,-1) + 3*data.move(index,0,0,1) - 6*data(index) + 3*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_central_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(index,3,1,0) + 4*data.move(index,2,1,0) - 5*data.move(index,1,1,0) + 2*data.move(index,0,1,0) + data.move(index,3,-1,0)
                 - 4*data.move(index,2,-1,0) + 5*data.move(index,1,-1,0) - 2*data.move(index,0,-1,0)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(index,1,3,0) + 4*data.move(index,1,2,0) - 5*data.move(index,1,1,0) + 2*data.move(index,1,0,0) + data.move(index,-1,3,0)
                - 4*data.move(index,-1,2,0) + 5*data.move(index,-1,1,0) - 2*data.move(index,-1,0,0)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(index,3,0,1) + 4*data.move(index,2,0,1) - 5*data.move(index,1,0,1) + 2*data.move(index,0,0,1) + data.move(index,3,0,-1)
                - 4*data.move(index,2,0,-1) + 5*data.move(index,1,0,-1) - 2*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(index,1,0,3) + 4*data.move(index,1,0,2) - 5*data.move(index,1,0,1) + 2*data.move(index,1,0,0) + data.move(index,-1,0,3)
                - 4*data.move(index,-1,0,2) + 5*data.move(index,-1,0,1) - 2*data.move(index,-1,0,0)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(index,0,3,1) + 4*data.move(index,0,2,1) - 5*data.move(index,0,1,1) + 2*data.move(index,0,0,1) + data.move(index,0,3,-1)
                - 4*data.move(index,0,2,-1) + 5*data.move(index,0,1,-1) - 2*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(index,0,1,3) + 4*data.move(index,0,1,2) - 5*data.move(index,0,1,1) + 2*data.move(index,0,1,0) + data.move(index,0,-1,3)
                - 4*data.move(index,0,-1,2) + 5*data.move(index,0,-1,1) - 2*data.move(index,0,-1,0)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto forward_backward_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(index,3,-2,0) + 4*data.move(index,2,-2,0) - 5*data.move(index,1,-2,0) + 2*data.move(index,0,-2,0) + 4*data.move(index,3,-1,0)
         - 16*data.move(index,2,-1,0) + 20*data.move(index,1,-1,0) - 8*data.move(index,0,-1,0) - 3*data.move(index,3,0,0) + 12*data.move(index,2,0,0)
         - 15*data.move(index,1,0,0) + 6*data(index)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(index,-2,3,0) + 4*data.move(index,-2,2,0) - 5*data.move(index,-2,1,0) + 2*data.move(index,-2,0,0) + 4*data.move(index,-1,3,0)
                - 16*data.move(index,-1,2,0) + 20*data.move(index,-1,1,0) - 8*data.move(index,-1,0,0) - 3*data.move(index,0,3,0) + 12*data.move(index,0,2,0)
                - 15*data.move(index,0,1,0) + 6*data(index)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(index,3,0,-2) + 4*data.move(index,2,0,-2) - 5*data.move(index,1,0,-2) + 2*data.move(index,0,0,-2) + 4*data.move(index,3,0,-1)
                - 16*data.move(index,2,0,-1) + 20*data.move(index,1,0,-1) - 8*data.move(index,0,0,-1) - 3*data.move(index,3,0,0) + 12*data.move(index,2,0,0)
                - 15*data.move(index,1,0,0) + 6*data(index)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(index,-2,0,3) + 4*data.move(index,-2,0,2) - 5*data.move(index,-2,0,1) + 2*data.move(index,-2,0,0) + 4*data.move(index,-1,0,3)
                - 16*data.move(index,-1,0,2) + 20*data.move(index,-1,0,1) - 8*data.move(index,-1,0,0) - 3*data.move(index,0,0,3) + 12*data.move(index,0,0,2)
                - 15*data.move(index,0,0,1) + 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(index,0,3,-2) + 4*data.move(index,0,2,-2) - 5*data.move(index,0,1,-2) + 2*data.move(index,0,0,-2) + 4*data.move(index,0,3,-1)
                - 16*data.move(index,0,2,-1) + 20*data.move(index,0,1,-1) - 8*data.move(index,0,0,-1) - 3*data.move(index,0,3,0) + 12*data.move(index,0,2,0)
                - 15*data.move(index,0,1,0) + 6*data(index)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(index,0,-2,3) + 4*data.move(index,0,-2,2) - 5*data.move(index,0,-2,1) + 2*data.move(index,0,-2,0) + 4*data.move(index,0,-1,3)
                - 16*data.move(index,0,-1,2) + 20*data.move(index,0,-1,1) - 8*data.move(index,0,-1,0) - 3*data.move(index,0,0,3) + 12*data.move(index,0,0,2)
                - 15*data.move(index,0,0,1) + 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_central_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(index,-3,1,0) + 4*data.move(index,-2,1,0) - 5*data.move(index,-1,1,0) + 2*data.move(index,0,1,0) + data.move(index,-3,-1,0)
                - 4*data.move(index,-2,-1,0) + 5*data.move(index,-1,-1,0) - 2*data.move(index,0,-1,0)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(index,1,-3,0) + 4*data.move(index,1,-2,0) - 5*data.move(index,1,-1,0) + 2*data.move(index,1,0,0) + data.move(index,-1,-3,0)
                - 4*data.move(index,-1,-2,0) + 5*data.move(index,-1,-1,0) - 2*data.move(index,-1,0,0)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(index,-3,0,1) + 4*data.move(index,-2,0,1) - 5*data.move(index,-1,0,1) + 2*data.move(index,0,0,1) + data.move(index,-3,0,-1)
                - 4*data.move(index,-2,0,-1) + 5*data.move(index,-1,0,-1) - 2*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(index,1,0,-3) + 4*data.move(index,1,0,-2) - 5*data.move(index,1,0,-1) + 2*data.move(index,1,0,0) + data.move(index,-1,0,-3)
                - 4*data.move(index,-1,0,-2) + 5*data.move(index,-1,0,-1) - 2*data.move(index,-1,0,0)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(index,0,-3,1) + 4*data.move(index,0,-2,1) - 5*data.move(index,0,-1,1) + 2*data.move(index,0,0,1) + data.move(index,0,-3,-1)
                - 4*data.move(index,0,-2,-1) + 5*data.move(index,0,-1,-1) - 2*data.move(index,0,0,-1)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(index,0,1,-3) + 4*data.move(index,0,1,-2) - 5*data.move(index,0,1,-1) + 2*data.move(index,0,1,0) + data.move(index,0,-1,-3)
                - 4*data.move(index,0,-1,-2) + 5*data.move(index,0,-1,-1) - 2*data.move(index,0,-1,0)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
#ifdef ALWAYS_INLINE_DERIVS
__attribute__((always_inline))
#endif
inline constexpr auto backward_forward_difference_3rd_mixed(const T& data, const size_t index)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(index,-3,2,0) - 4*data.move(index,-2,2,0) + 5*data.move(index,-1,2,0) - 2*data.move(index,0,2,0) - 4*data.move(index,-3,1,0)
                + 16*data.move(index,-2,1,0) - 20*data.move(index,-1,1,0) + 8*data.move(index,0,1,0) + 3*data.move(index,-3,0,0) - 12*data.move(index,-2,0,0)
                + 15*data.move(index,-1,0,0) - 6*data(index)) / (2*data.dy(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(index,2,-3,0) - 4*data.move(index,2,-2,0) + 5*data.move(index,2,-1,0) - 2*data.move(index,2,0,0) - 4*data.move(index,1,-3,0)
                + 16*data.move(index,1,-2,0) - 20*data.move(index,1,-1,0) + 8*data.move(index,1,0,0) + 3*data.move(index,0,-3,0) - 12*data.move(index,0,-2,0)
                + 15*data.move(index,0,-1,0) - 6*data(index)) / (2*data.dy(index)*data.dy(index)*data.dx(index));
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(index,-3,0,2) - 4*data.move(index,-2,0,2) + 5*data.move(index,-1,0,2) - 2*data.move(index,0,0,2) - 4*data.move(index,-3,0,1)
                + 16*data.move(index,-2,0,1) - 20*data.move(index,-1,0,1) + 8*data.move(index,0,0,1) + 3*data.move(index,-3,0,0) - 12*data.move(index,-2,0,0)
                + 15*data.move(index,-1,0,0) - 6*data(index)) / (2*data.dz(index)*data.dx(index)*data.dx(index));
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(index,2,0,-3) - 4*data.move(index,2,0,-2) + 5*data.move(index,2,0,-1) - 2*data.move(index,2,0,0) - 4*data.move(index,1,0,-3)
                + 16*data.move(index,1,0,-2) - 20*data.move(index,1,0,-1) + 8*data.move(index,1,0,0) + 3*data.move(index,0,0,-3) - 12*data.move(index,0,0,-2)
                + 15*data.move(index,0,0,-1) - 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dx(index));
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(index,0,-3,2) - 4*data.move(index,0,-2,2) + 5*data.move(index,0,-1,2) - 2*data.move(index,0,0,2) - 4*data.move(index,0,-3,1)
                + 16*data.move(index,0,-2,1) - 20*data.move(index,0,-1,1) + 8*data.move(index,0,0,1) + 3*data.move(index,0,-3,0) - 12*data.move(index,0,-2,0)
                + 15*data.move(index,0,-1,0) - 6*data(index)) / (2*data.dz(index)*data.dy(index)*data.dy(index));
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(index,0,2,-3) - 4*data.move(index,0,2,-2) + 5*data.move(index,0,2,-1) - 2*data.move(index,0,2,0) - 4*data.move(index,0,1,-3)
                + 16*data.move(index,0,1,-2) - 20*data.move(index,0,1,-1) + 8*data.move(index,0,1,0) + 3*data.move(index,0,0,-3) - 12*data.move(index,0,0,-2)
                + 15*data.move(index,0,0,-1) - 6*data(index)) / (2*data.dz(index)*data.dz(index)*data.dy(index));
    }
}



#endif //CODE_FINITE_DIFFERENCE_HPP
