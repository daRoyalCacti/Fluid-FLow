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
inline constexpr auto central_difference_1st(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data.move(x,y,z,1,0,0) - data.move(x,y,z,-1,0,0) ) / (2*dx);
    } else if constexpr(axis == 1) {
        return ( data.move(x,y,z,0,1,0) - data.move(x,y,z,0,-1,0) ) / (2*dy);
    } else if constexpr(axis == 2) {
        return ( data.move(x,y,z,0,0,1) - data.move(x,y,z,0,0,-1) ) / (2*dz);
    }
}

template <unsigned axis, typename T>
inline constexpr auto forward_difference_1st(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( -data.move(x,y,z,2,0,0) + 4*data.move(x,y,z,1,0,0) - 3*data(x,y,z)) / (2*dx);
    } else if constexpr(axis == 1) {
        return ( -data.move(x,y,z,0,2,0) + 4*data.move(x,y,z,0,1,0) - 3*data(x,y,z)) / (2*dy);
    } else if constexpr(axis == 2) {
        return ( -data.move(x,y,z,0,0,2) + 4*data.move(x,y,z,0,0,1) - 3*data(x,y,z)) / (2*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto backward_difference_1st(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( 3*data(x,y,z) - 4*data.move(x,y,z,-1,0,0) + data.move(x,y,z,-2,0,0) ) / (2*dx);
    } else if constexpr(axis == 1) {
        return ( 3*data(x,y,z) - 4*data.move(x,y,z,0,-1,0) + data.move(x,y,z,0,-2,0) ) / (2*dy);
    } else if constexpr(axis == 2) {
        return ( 3*data(x,y,z) - 4*data.move(x,y,z,0,0,-1) + data.move(x,y,z,0,0,-2) ) / (2*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto central_difference_2nd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data.move(x,y,z,1,0,0) - 2* data(x,y,z) + data.move(x,y,z,-1,0,0) ) / (dx*dx);
    } else if constexpr(axis == 1) {
        return  ( data.move(x,y,z,0,1,0) - 2* data(x,y,z) + data.move(x,y,z,0,-1,0) ) / (dy*dy);
    } else if constexpr(axis == 2) {
        return  ( data.move(x,y,z,0,0,1) - 2* data(x,y,z) + data.move(x,y,z,0,0,-1) ) / (dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto forward_difference_2nd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (-data.move(x,y,z,3,0,0)+4*data.move(x,y,z,2,0,0) - 5*data.move(x,y,z,1,0,0)+2*data(x,y,z)) / (dx*dx);
    } else if constexpr(axis == 1) {
        return (-data.move(x,y,z,0,3,0)+4*data.move(x,y,z,0,2,0) - 5*data.move(x,y,z,0,1,0)+2*data(x,y,z)) / (dy*dy);
    } else if constexpr(axis == 2) {
        return (-data.move(x,y,z,0,0,3)+4*data.move(x,y,z,0,0,2) - 5*data.move(x,y,z,0,0,1)+2*data(x,y,z)) / (dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto backward_difference_2nd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( 2*data(x,y,z) -5*data.move(x,y,z,-1,0,0)+4*data.move(x,y,z,-2,0,0) - data.move(x,y,z,-3,0,0) ) / (dx*dx);
    } else if constexpr(axis == 1) {
        return ( 2*data(x,y,z) -5*data.move(x,y,z,0,-1,0)+4*data.move(x,y,z,0,-2,0) - data.move(x,y,z,0,-3,0) ) / (dy*dy);
    } else if constexpr(axis == 2) {
        return ( 2*data(x,y,z) -5*data.move(x,y,z,0,0,-1)+4*data.move(x,y,z,0,0,-2) - data.move(x,y,z,0,0,-3) ) / (dz*dz);
    }
}

template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return central_difference_2nd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data.move(x,y,z,1,1,0) - data.move(x,y,z,-1,1,0) - data.move(x,y,z,1,-1,0) + data(x-1,y-1,z))  / (4*dx*dy);
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data.move(x,y,z,1,0,1) - data.move(x,y,z,-1,0,1) - data.move(x,y,z,1,0,-1) + data.move(x,y,z,-1,0,-1))  / (4*dx*dz);
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data.move(x,y,z,0,1,1) - data.move(x,y,z,0,-1,1) - data.move(x,y,z,0,1,-1) + data.move(x,y,z,0,-1,-1))  / (4*dy*dz);
    }
}

template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return forward_difference_2nd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data.move(x,y,z,2,2,0) - 4*data.move(x,y,z,1,2,0) + 3*data.move(x,y,z,0,2,0) - 4*data.move(x,y,z,2,1,0)+16*data.move(x,y,z,1,1,0)
                    -12*data.move(x,y,z,0,1,0) + 3*data.move(x,y,z,2,0,0) - 12*data.move(x,y,z,1,0,0) + 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data.move(x,y,z,2,0,2) - 4*data.move(x,y,z,1,0,2) + 3*data.move(x,y,z,0,0,2) - 4*data.move(x,y,z,2,0,1)+16*data.move(x,y,z,1,0,1)
                 -12*data.move(x,y,z,0,0,1) + 3*data.move(x,y,z,2,0,0) - 12*data.move(x,y,z,1,0,0) + 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data.move(x,y,z,0,2,2) - 4*data.move(x,y,z,0,1,2) + 3*data.move(x,y,z,0,0,2) - 4*data.move(x,y,z,0,2,1)+16*data.move(x,y,z,0,1,1)
                 -12*data.move(x,y,z,0,0,1) + 3*data.move(x,y,z,0,2,0) - 12*data.move(x,y,z,0,1,0) + 9*data(x,y,z)) / (4*dy*dz);
    }
}


template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return backward_difference_2nd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data.move(x,y,z,-2,-2,0) - 4*data.move(x,y,z,-1,-2,0) + 3*data.move(x,y,z,0,-2,0) - 4*data.move(x,y,z,-2,-1,0)+16*data.move(x,y,z,-1,-1,0)
                 -12*data.move(x,y,z,0,-1,0) + 3*data.move(x,y,z,-2,0,0) - 12*data.move(x,y,z,-1,0,0) + 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data.move(x,y,z,-2,0,-2) - 4*data.move(x,y,z,-1,0,-2) + 3*data.move(x,y,z,0,0,-2) - 4*data.move(x,y,z,-2,0,-1)+16*data.move(x,y,z,-1,0,-1)
                 -12*data.move(x,y,z,0,0,-1) + 3*data.move(x,y,z,-2,0,0) - 12*data.move(x,y,z,-1,0,0) + 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data.move(x,y,z,0,-2,-2) - 4*data.move(x,y,z,0,-1,-2) + 3*data.move(x,y,z,0,0,-2) - 4*data.move(x,y,z,0,-2,-1)+16*data.move(x,y,z,0,-1,-1)
                 -12*data.move(x,y,z,0,0,-1) + 3*data.move(x,y,z,0,-2,0) - 12*data.move(x,y,z,0,-1,0) + 9*data(x,y,z)) / (4*dy*dz);
    }
}



template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_forward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return (-data.move(x,y,z,2,1,0) + 4*data.move(x,y,z,1,1,0) - 3*data.move(x,y,z,0,1,0) + data.move(x,y,z,2,-1,0) - 4*data.move(x,y,z,1,-1,0)+3*data.move(x,y,z,0,-1,0)) / (4*dx*dy);
    } else if constexpr(axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return (-data.move(x,y,z,1,2,0) + 4*data.move(x,y,z,1,1,0) - 3*data.move(x,y,z,1,0,0) + data.move(x,y,z,-1,2,0) - 4*data.move(x,y,z,-1,1,0)+3*data.move(x,y,z,-1,0,0)) / (4*dx*dy);
    } else if constexpr( axis1 == 2 && axis2 == 0) {//d^2/dxdz
        return (-data.move(x,y,z,2,0,1) + 4*data.move(x,y,z,1,0,1) - 3*data.move(x,y,z,0,0,1) + data.move(x,y,z,2,0,-1) - 4*data.move(x,y,z,1,0,-1)+3*data.move(x,y,z,0,0,-1)) / (4*dx*dz);
    } else if constexpr(axis1 == 0 && axis2 == 2 ) {  //d^2/dxdz
        return (-data.move(x,y,z,1,0,2) + 4*data.move(x,y,z,1,0,1) - 3*data.move(x,y,z,1,0,0) + data.move(x,y,z,-1,0,2) - 4*data.move(x,y,z,-1,0,1)+3*data.move(x,y,z,-1,0,0)) / (4*dx*dz);
    } else if constexpr( axis1 == 2 && axis2 == 1) {//d^2/dydz
        return (-data.move(x,y,z,0,2,1) + 4*data.move(x,y,z,0,1,1) - 3*data.move(x,y,z,0,0,1) + data.move(x,y,z,0,2,-1) - 4*data.move(x,y,z,0,1,-1)+3*data.move(x,y,z,0,0,-1)) / (4*dy*dz);
    } else if constexpr(axis1 == 1 && axis2 == 2) { //d^2/dydz
        return (-data.move(x,y,z,0,1,2) + 4*data.move(x,y,z,0,1,1) - 3*data.move(x,y,z,0,1,0) + data.move(x,y,z,0,-1,2) - 4*data.move(x,y,z,0,-1,1)+3*data.move(x,y,z,0,-1,0)) / (4*dy*dz);
    }
}


template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_backward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return (data.move(x,y,z,-2,1,0) - 4*data.move(x,y,z,-1,1,0) + 3*data.move(x,y,z,0,1,0) - data.move(x,y,z,-2,-1,0) + 4*data.move(x,y,z,-1,-1,0)-3*data.move(x,y,z,0,-1,0)) / (4*dx*dy);
    } else if constexpr(axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return (data.move(x,y,z,1,-2,0) - 4*data.move(x,y,z,1,-1,0) + 3*data.move(x,y,z,1,0,0) - data.move(x,y,z,-1,-2,0) + 4*data.move(x,y,z,-1,-1,0)-3*data.move(x,y,z,-1,0,0)) / (4*dx*dy);
    } else if constexpr( axis1 == 2 && axis2 == 0) {//d^2/dxdz
        return (data.move(x,y,z,-2,0,1) - 4*data.move(x,y,z,-1,0,1) + 3*data.move(x,y,z,0,0,1) - data.move(x,y,z,-2,0,-1) + 4*data.move(x,y,z,-1,0,-1)-3*data.move(x,y,z,0,0,-1)) / (4*dx*dz);
    } else if constexpr(axis1 == 0 && axis2 == 2 ) {  //d^2/dxdz
        return (data.move(x,y,z,1,0,-2) - 4*data.move(x,y,z,1,0,-1) + 3*data.move(x,y,z,1,0,0) - data.move(x,y,z,-1,0,-2) + 4*data.move(x,y,z,-1,0,-1)-3*data.move(x,y,z,-1,0,0)) / (4*dx*dz);
    } else if constexpr( axis1 == 2 && axis2 == 1) {//d^2/dydz
        return (data.move(x,y,z,0,-2,1) - 4*data.move(x,y,z,0,-1,1) + 3*data.move(x,y,z,0,0,1) - data.move(x,y,z,0,-2,-1) + 4*data.move(x,y,z,0,-1,-1)-3*data.move(x,y,z,0,0,-1)) / (4*dy*dz);
    } else if constexpr(axis1 == 1 && axis2 == 2) { //d^2/dydz
        return (data.move(x,y,z,0,1,-2) - 4*data.move(x,y,z,0,1,-1) + 3*data.move(x,y,z,0,1,0) - data.move(x,y,z,0,-1,-2) + 4*data.move(x,y,z,0,-1,-1)-3*data.move(x,y,z,0,-1,0)) / (4*dy*dz);
    }
}


template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_backward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return ( -data.move(x,y,z,2,-2,0) + 4*data.move(x,y,z,1,-2,0) - 3*data.move(x,y,z,0,-2,0) + 4*data.move(x,y,z,2,-1,0)-16*data.move(x,y,z,1,-1,0)
                 +12*data.move(x,y,z,0,-1,0) - 3*data.move(x,y,z,2,0,0) + 12*data.move(x,y,z,1,0,0) - 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return ( -data.move(x,y,z,-2,2,0) + 4*data.move(x,y,z,-1,2,0) - 3*data.move(x,y,z,0,2,0) + 4*data.move(x,y,z,-2,1,0)-16*data.move(x,y,z,-1,1,0)
                 +12*data.move(x,y,z,0,1,0) - 3*data.move(x,y,z,-2,0,0) + 12*data.move(x,y,z,-1,0,0) - 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^2/dxdz
        return ( -data.move(x,y,z,2,0,-2) + 4*data.move(x,y,z,1,0,-2) - 3*data.move(x,y,z,0,0,-2) + 4*data.move(x,y,z,2,0,-1)-16*data.move(x,y,z,1,0,-1)
                 +12*data.move(x,y,z,0,0,-1) - 3*data.move(x,y,z,2,0,0) + 12*data.move(x,y,z,1,0,0) - 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^2/dxdz
        return ( -data.move(x,y,z,-2,0,2) + 4*data.move(x,y,z,-1,0,2) - 3*data.move(x,y,z,0,0,2) + 4*data.move(x,y,z,-2,0,1)-16*data.move(x,y,z,-1,0,1)
                 +12*data.move(x,y,z,0,0,1) - 3*data.move(x,y,z,-2,0,0) + 12*data.move(x,y,z,-1,0,0) - 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^2/dydz
        return ( -data.move(x,y,z,0,2,-2) + 4*data.move(x,y,z,0,1,-2) - 3*data.move(x,y,z,0,0,-2) + 4*data.move(x,y,z,0,2,-1)-16*data.move(x,y,z,0,1,-1)
                 +12*data.move(x,y,z,0,0,-1) - 3*data.move(x,y,z,0,2,0) + 12*data.move(x,y,z,0,1,0) - 9*data(x,y,z)) / (4*dy*dz);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^2/dydz
        return ( -data.move(x,y,z,0,-2,2) + 4*data.move(x,y,z,0,-1,2) - 3*data.move(x,y,z,0,0,2) + 4*data.move(x,y,z,0,-2,1)-16*data.move(x,y,z,0,-1,1)
                 +12*data.move(x,y,z,0,0,1) - 3*data.move(x,y,z,0,-2,0) + 12*data.move(x,y,z,0,-1,0) - 9*data(x,y,z)) / (4*dy*dz);
    }
}



template <unsigned axis, typename T>
inline constexpr auto central_difference_3rd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data.move(x,y,z,2,0,0) - 2*data.move(x,y,z,1,0,0) + 2*data.move(x,y,z,-1,0,0) - data.move(x,y,z,-2,0,0)) / (2*dx*dx*dx);
    } else if constexpr(axis == 1) {
        return ( data.move(x,y,z,0,2,0) - 2*data.move(x,y,z,0,1,0) + 2*data.move(x,y,z,0,-1,0) - data.move(x,y,z,0,-2,0)) / (2*dy*dy*dy);
    } else if constexpr(axis == 2) {
        return ( data.move(x,y,z,0,0,2) - 2*data.move(x,y,z,0,0,1) + 2*data.move(x,y,z,0,0,-1) - data.move(x,y,z,0,0,-2)) / (2*dz*dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto forward_difference_3rd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (-3*data.move(x,y,z,4,0,0)+14*data.move(x,y,z,3,0,0) - 24*data.move(x,y,z,2,0,0)+ 18*data.move(x,y,z,1,0,0) - 5*data(x,y,z)) / (2*dx*dx*dx);
    } else if constexpr(axis == 1) {
        return (-3*data.move(x,y,z,0,4,0)+14*data.move(x,y,z,0,3,0) - 24*data.move(x,y,z,0,2,0)+ 18*data.move(x,y,z,0,1,0) - 5*data(x,y,z)) / (2*dy*dy*dy);
    } else if constexpr(axis == 2) {
        return (-3*data.move(x,y,z,0,0,4)+14*data.move(x,y,z,0,0,3) - 24*data.move(x,y,z,0,0,2)+ 18*data.move(x,y,z,0,0,1) - 5*data(x,y,z)) / (2*dz*dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto backward_difference_3rd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (5*data(x,y,z) - 18*data.move(x,y,z,-1,0,0) + 24*data.move(x,y,z,-2,0,0) - 14*data.move(x,y,z,-3,0,0) + 3*data.move(x,y,z,-4,0,0)) / (2*dx*dx*dx);
    } else if constexpr(axis == 1) {
        return (5*data(x,y,z) - 18*data.move(x,y,z,0,-1,0) + 24*data.move(x,y,z,0,-2,0) - 14*data.move(x,y,z,0,-3,0) + 3*data.move(x,y,z,0,-4,0)) / (2*dy*dy*dy);
    } else if constexpr(axis == 2) {
        return (5*data(x,y,z) - 18*data.move(x,y,z,0,0,-1) + 24*data.move(x,y,z,0,0,-2) - 14*data.move(x,y,z,0,0,-3) + 3*data.move(x,y,z,0,0,-4)) / (2*dz*dz*dz);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return central_difference_3rd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(x,y,z,1,1,0) - 2*data.move(x,y,z,0,1,0) + data.move(x,y,z,-1,1,0) - data.move(x,y,z,1,-1,0)+2*data.move(x,y,z,0,-1,0)-data.move(x,y,z,-1,-1,0)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(x,y,z,1,1,0) - 2*data.move(x,y,z,1,0,0) + data.move(x,y,z,1,-1,0) - data.move(x,y,z,-1,1,0)+2*data.move(x,y,z,-1,0,0)-data.move(x,y,z,-1,-1,0)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(x,y,z,1,0,1) - 2*data.move(x,y,z,0,0,1) + data.move(x,y,z,-1,0,1) - data.move(x,y,z,1,0,-1)+2*data.move(x,y,z,0,0,-1)-data.move(x,y,z,-1,0,-1)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(x,y,z,1,0,1) - 2*data.move(x,y,z,1,0,0) + data.move(x,y,z,1,0,-1) - data.move(x,y,z,-1,0,1)+2*data.move(x,y,z,-1,0,0)-data.move(x,y,z,-1,0,-1)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(x,y,z,0,1,1) - 2*data.move(x,y,z,0,0,1) + data.move(x,y,z,0,-1,1) - data.move(x,y,z,0,1,-1)+2*data.move(x,y,z,0,0,-1)-data.move(x,y,z,0,-1,-1)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(x,y,z,0,1,1) - 2*data.move(x,y,z,0,1,0) + data.move(x,y,z,0,1,-1) - data.move(x,y,z,0,-1,1)+2*data.move(x,y,z,0,-1,0)-data.move(x,y,z,0,-1,-1)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return forward_difference_3rd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(x,y,z,3,2,0) - 4*data.move(x,y,z,2,2,0) + 5*data.move(x,y,z,1,2,0) - 2*data.move(x,y,z,0,2,0) - 4*data.move(x,y,z,3,1,0)
                 + 16*data.move(x,y,z,2,1,0) - 20*data.move(x,y,z,1,1,0) + 8*data.move(x,y,z,0,1,0) + 3*data.move(x,y,z,3,0,0) - 12*data.move(x,y,z,2,0,0)
                 + 15*data.move(x,y,z,1,0,0) - 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(x,y,z,2,3,0) - 4*data.move(x,y,z,2,2,0) + 5*data.move(x,y,z,2,1,0) - 2*data.move(x,y,z,2,0,0) - 4*data.move(x,y,z,1,3,0)
                + 16*data.move(x,y,z,1,2,0) - 20*data.move(x,y,z,1,1,0) + 8*data.move(x,y,z,1,0,0) + 3*data.move(x,y,z,0,3,0) - 12*data.move(x,y,z,0,2,0)
                + 15*data.move(x,y,z,0,1,0) - 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(x,y,z,3,0,2) - 4*data.move(x,y,z,2,0,2) + 5*data.move(x,y,z,1,0,2) - 2*data.move(x,y,z,0,0,2) - 4*data.move(x,y,z,3,0,1)
                + 16*data.move(x,y,z,2,0,1) - 20*data.move(x,y,z,1,0,1) + 8*data.move(x,y,z,0,0,1) + 3*data.move(x,y,z,3,0,0) - 12*data.move(x,y,z,2,0,0)
                + 15*data.move(x,y,z,1,0,0) - 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(x,y,z,2,0,3) - 4*data.move(x,y,z,2,0,2) + 5*data.move(x,y,z,2,0,1) - 2*data.move(x,y,z,2,0,0) - 4*data.move(x,y,z,1,0,3)
                + 16*data.move(x,y,z,1,0,2) - 20*data.move(x,y,z,1,0,1) + 8*data.move(x,y,z,1,0,0) + 3*data.move(x,y,z,0,0,3) - 12*data.move(x,y,z,0,0,2)
                + 15*data.move(x,y,z,0,0,1) - 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(x,y,z,0,3,2) - 4*data.move(x,y,z,0,2,2) + 5*data.move(x,y,z,0,1,2) - 2*data.move(x,y,z,0,0,2) - 4*data.move(x,y,z,0,3,1)
                + 16*data.move(x,y,z,0,2,1) - 20*data.move(x,y,z,0,1,1) + 8*data.move(x,y,z,0,0,1) + 3*data.move(x,y,z,0,3,0) - 12*data.move(x,y,z,0,2,0)
                + 15*data.move(x,y,z,0,1,0) - 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(x,y,z,0,2,3) - 4*data.move(x,y,z,0,2,2) + 5*data.move(x,y,z,0,2,1) - 2*data.move(x,y,z,0,2,0) - 4*data.move(x,y,z,0,1,3)
                + 16*data.move(x,y,z,0,1,2) - 20*data.move(x,y,z,0,1,1) + 8*data.move(x,y,z,0,1,0) + 3*data.move(x,y,z,0,0,3) - 12*data.move(x,y,z,0,0,2)
                + 15*data.move(x,y,z,0,0,1) - 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return backward_difference_3rd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(x,y,z,-3,-2,0) + 4*data.move(x,y,z,-2,-2,0) - 5*data.move(x,y,z,-1,-2,0) + 2*data.move(x,y,z,0,-2,0) + 4*data.move(x,y,z,-3,-1,0)
                - 16*data.move(x,y,z,-2,-1,0) + 20*data.move(x,y,z,-1,-1,0) - 8*data.move(x,y,z,0,-1,0) - 3*data.move(x,y,z,-3,0,0) + 12*data.move(x,y,z,-2,0,0)
                - 15*data.move(x,y,z,-1,0,0) + 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(x,y,z,-2,-3,0) + 4*data.move(x,y,z,-2,-2,0) - 5*data.move(x,y,z,-2,-1,0) + 2*data.move(x,y,z,-2,0,0) + 4*data.move(x,y,z,-1,-3,0)
                - 16*data.move(x,y,z,-1,-2,0) + 20*data.move(x,y,z,-1,-1,0) - 8*data.move(x,y,z,-1,0,0) - 3*data.move(x,y,z,0,-3,0) + 12*data.move(x,y,z,0,-2,0)
                - 15*data.move(x,y,z,0,-1,0) + 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(x,y,z,-3,0,-2) + 4*data.move(x,y,z,-2,0,-2) - 5*data.move(x,y,z,-1,0,-2) + 2*data.move(x,y,z,0,0,-2) + 4*data.move(x,y,z,-3,0,-1)
                - 16*data.move(x,y,z,-2,0,-1) + 20*data.move(x,y,z,-1,0,-1) - 8*data.move(x,y,z,0,0,-1) - 3*data.move(x,y,z,-3,0,0) + 12*data.move(x,y,z,-2,0,0)
                - 15*data.move(x,y,z,-1,0,0) + 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(x,y,z,-2,0,-3) + 4*data.move(x,y,z,-2,0,-2) - 5*data.move(x,y,z,-2,0,-1) + 2*data.move(x,y,z,-2,0,0) + 4*data.move(x,y,z,-1,0,-3)
                - 16*data.move(x,y,z,-1,0,-2) + 20*data.move(x,y,z,-1,0,-1) - 8*data.move(x,y,z,-1,0,0) - 3*data.move(x,y,z,0,0,-3) + 12*data.move(x,y,z,0,0,-2)
                - 15*data.move(x,y,z,0,0,-1) + 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(x,y,z,0,-3,-2) + 4*data.move(x,y,z,0,-2,-2) - 5*data.move(x,y,z,0,-1,-2) + 2*data.move(x,y,z,0,0,-2) + 4*data.move(x,y,z,0,-3,-1)
                - 16*data.move(x,y,z,0,-2,-1) + 20*data.move(x,y,z,0,-1,-1) - 8*data.move(x,y,z,0,0,-1) - 3*data.move(x,y,z,0,-3,0) + 12*data.move(x,y,z,0,-2,0)
                - 15*data.move(x,y,z,0,-1,0) + 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(x,y,z,0,-2,-3) + 4*data.move(x,y,z,0,-2,-2) - 5*data.move(x,y,z,0,-2,-1) + 2*data.move(x,y,z,0,-2,0) + 4*data.move(x,y,z,0,-1,-3)
                - 16*data.move(x,y,z,0,-1,-2) + 20*data.move(x,y,z,0,-1,-1) - 8*data.move(x,y,z,0,-1,0) - 3*data.move(x,y,z,0,0,-3) + 12*data.move(x,y,z,0,0,-2)
                - 15*data.move(x,y,z,0,0,-1) + 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_forward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(x,y,z,1,2,0) + 2*data.move(x,y,z,0,2,0) - data.move(x,y,z,-1,2,0) + 4*data.move(x,y,z,1,1,0) - 8*data.move(x,y,z,0,1,0)
                 + 4*data.move(x,y,z,-1,1,0) - 3*data.move(x,y,z,1,0,0) + 6*data(x,y,z) - 3*data.move(x,y,z,-1,0,0)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(x,y,z,2,1,0) + 2*data.move(x,y,z,2,0,0) - data.move(x,y,z,2,-1,0) + 4*data.move(x,y,z,1,1,0) - 8*data.move(x,y,z,1,0,0)
                + 4*data.move(x,y,z,1,-1,0) - 3*data.move(x,y,z,0,1,0) + 6*data(x,y,z) - 3*data.move(x,y,z,0,-1,0)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(x,y,z,1,0,2) + 2*data.move(x,y,z,0,0,2) - data.move(x,y,z,-1,0,2) + 4*data.move(x,y,z,1,0,1) - 8*data.move(x,y,z,0,0,1)
                + 4*data.move(x,y,z,-1,0,1) - 3*data.move(x,y,z,1,0,0) + 6*data(x,y,z) - 3*data.move(x,y,z,-1,0,0)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(x,y,z,2,0,1) + 2*data.move(x,y,z,2,0,0) - data.move(x,y,z,2,0,-1) + 4*data.move(x,y,z,1,0,1) - 8*data.move(x,y,z,1,0,0)
                + 4*data.move(x,y,z,1,0,-1) - 3*data.move(x,y,z,0,0,1) + 6*data(x,y,z) - 3*data.move(x,y,z,0,0,-1)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(x,y,z,0,1,2) + 2*data.move(x,y,z,0,0,2) - data.move(x,y,z,0,-1,2) + 4*data.move(x,y,z,0,1,1) - 8*data.move(x,y,z,0,0,1)
                + 4*data.move(x,y,z,0,-1,1) - 3*data.move(x,y,z,0,1,0) + 6*data(x,y,z) - 3*data.move(x,y,z,0,-1,0)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(x,y,z,0,2,1) + 2*data.move(x,y,z,0,2,0) - data.move(x,y,z,0,2,-1) + 4*data.move(x,y,z,0,1,1) - 8*data.move(x,y,z,0,1,0)
                + 4*data.move(x,y,z,0,1,-1) - 3*data.move(x,y,z,0,0,1) + 6*data(x,y,z) - 3*data.move(x,y,z,0,0,-1)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_backward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(x,y,z,1,-2,0) - 2*data.move(x,y,z,0,-2,0) + data.move(x,y,z,-1,-2,0) - 4*data.move(x,y,z,1,-1,0) + 8*data.move(x,y,z,0,-1,0)
                - 4*data.move(x,y,z,-1,-1,0) + 3*data.move(x,y,z,1,0,0) - 6*data(x,y,z) + 3*data.move(x,y,z,-1,0,0)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(x,y,z,-2,1,0) - 2*data.move(x,y,z,-2,0,0) + data.move(x,y,z,-2,-1,0) - 4*data.move(x,y,z,-1,1,0) + 8*data.move(x,y,z,-1,0,0)
                - 4*data.move(x,y,z,-1,-1,0) + 3*data.move(x,y,z,0,1,0) - 6*data(x,y,z) + 3*data.move(x,y,z,0,-1,0)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(x,y,z,1,0,-2) - 2*data.move(x,y,z,0,0,-2) + data.move(x,y,z,-1,0,-2) - 4*data.move(x,y,z,1,0,-1) + 8*data.move(x,y,z,0,0,-1)
                - 4*data.move(x,y,z,-1,0,-1) + 3*data.move(x,y,z,1,0,0) - 6*data(x,y,z) + 3*data.move(x,y,z,-1,0,0)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(x,y,z,-2,0,1) - 2*data.move(x,y,z,-2,0,0) + data.move(x,y,z,-2,0,-1) - 4*data.move(x,y,z,-1,0,1) + 8*data.move(x,y,z,-1,0,0)
                - 4*data.move(x,y,z,-1,0,-1) + 3*data.move(x,y,z,0,0,1) - 6*data(x,y,z) + 3*data.move(x,y,z,0,0,-1)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(x,y,z,0,1,-2) - 2*data.move(x,y,z,0,0,-2) + data.move(x,y,z,0,-1,-2) - 4*data.move(x,y,z,0,1,-1) + 8*data.move(x,y,z,0,0,-1)
                - 4*data.move(x,y,z,0,-1,-1) + 3*data.move(x,y,z,0,1,0) - 6*data(x,y,z) + 3*data.move(x,y,z,0,-1,0)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(x,y,z,0,-2,1) - 2*data.move(x,y,z,0,-2,0) + data.move(x,y,z,0,-2,-1) - 4*data.move(x,y,z,0,-1,1) + 8*data.move(x,y,z,0,-1,0)
                - 4*data.move(x,y,z,0,-1,-1) + 3*data.move(x,y,z,0,0,1) - 6*data(x,y,z) + 3*data.move(x,y,z,0,0,-1)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_central_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(x,y,z,3,1,0) + 4*data.move(x,y,z,2,1,0) - 5*data.move(x,y,z,1,1,0) + 2*data.move(x,y,z,0,1,0) + data.move(x,y,z,3,-1,0)
                 - 4*data.move(x,y,z,2,-1,0) + 5*data.move(x,y,z,1,-1,0) - 2*data.move(x,y,z,0,-1,0)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(x,y,z,1,3,0) + 4*data.move(x,y,z,1,2,0) - 5*data.move(x,y,z,1,1,0) + 2*data.move(x,y,z,1,0,0) + data.move(x,y,z,-1,3,0)
                - 4*data.move(x,y,z,-1,2,0) + 5*data.move(x,y,z,-1,1,0) - 2*data.move(x,y,z,-1,0,0)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(x,y,z,3,0,1) + 4*data.move(x,y,z,2,0,1) - 5*data.move(x,y,z,1,0,1) + 2*data.move(x,y,z,0,0,1) + data.move(x,y,z,3,0,-1)
                - 4*data.move(x,y,z,2,0,-1) + 5*data.move(x,y,z,1,0,-1) - 2*data.move(x,y,z,0,0,-1)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(x,y,z,1,0,3) + 4*data.move(x,y,z,1,0,2) - 5*data.move(x,y,z,1,0,1) + 2*data.move(x,y,z,1,0,0) + data.move(x,y,z,-1,0,3)
                - 4*data.move(x,y,z,-1,0,2) + 5*data.move(x,y,z,-1,0,1) - 2*data.move(x,y,z,-1,0,0)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(x,y,z,0,3,1) + 4*data.move(x,y,z,0,2,1) - 5*data.move(x,y,z,0,1,1) + 2*data.move(x,y,z,0,0,1) + data.move(x,y,z,0,3,-1)
                - 4*data.move(x,y,z,0,2,-1) + 5*data.move(x,y,z,0,1,-1) - 2*data.move(x,y,z,0,0,-1)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(x,y,z,0,1,3) + 4*data.move(x,y,z,0,1,2) - 5*data.move(x,y,z,0,1,1) + 2*data.move(x,y,z,0,1,0) + data.move(x,y,z,0,-1,3)
                - 4*data.move(x,y,z,0,-1,2) + 5*data.move(x,y,z,0,-1,1) - 2*data.move(x,y,z,0,-1,0)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_backward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(x,y,z,3,-2,0) + 4*data.move(x,y,z,2,-2,0) - 5*data.move(x,y,z,1,-2,0) + 2*data.move(x,y,z,0,-2,0) + 4*data.move(x,y,z,3,-1,0)
         - 16*data.move(x,y,z,2,-1,0) + 20*data.move(x,y,z,1,-1,0) - 8*data.move(x,y,z,0,-1,0) - 3*data.move(x,y,z,3,0,0) + 12*data.move(x,y,z,2,0,0)
         - 15*data.move(x,y,z,1,0,0) + 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(x,y,z,-2,3,0) + 4*data.move(x,y,z,-2,2,0) - 5*data.move(x,y,z,-2,1,0) + 2*data.move(x,y,z,-2,0,0) + 4*data.move(x,y,z,-1,3,0)
                - 16*data.move(x,y,z,-1,2,0) + 20*data.move(x,y,z,-1,1,0) - 8*data.move(x,y,z,-1,0,0) - 3*data.move(x,y,z,0,3,0) + 12*data.move(x,y,z,0,2,0)
                - 15*data.move(x,y,z,0,1,0) + 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(x,y,z,3,0,-2) + 4*data.move(x,y,z,2,0,-2) - 5*data.move(x,y,z,1,0,-2) + 2*data.move(x,y,z,0,0,-2) + 4*data.move(x,y,z,3,0,-1)
                - 16*data.move(x,y,z,2,0,-1) + 20*data.move(x,y,z,1,0,-1) - 8*data.move(x,y,z,0,0,-1) - 3*data.move(x,y,z,3,0,0) + 12*data.move(x,y,z,2,0,0)
                - 15*data.move(x,y,z,1,0,0) + 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(x,y,z,-2,0,3) + 4*data.move(x,y,z,-2,0,2) - 5*data.move(x,y,z,-2,0,1) + 2*data.move(x,y,z,-2,0,0) + 4*data.move(x,y,z,-1,0,3)
                - 16*data.move(x,y,z,-1,0,2) + 20*data.move(x,y,z,-1,0,1) - 8*data.move(x,y,z,-1,0,0) - 3*data.move(x,y,z,0,0,3) + 12*data.move(x,y,z,0,0,2)
                - 15*data.move(x,y,z,0,0,1) + 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(x,y,z,0,3,-2) + 4*data.move(x,y,z,0,2,-2) - 5*data.move(x,y,z,0,1,-2) + 2*data.move(x,y,z,0,0,-2) + 4*data.move(x,y,z,0,3,-1)
                - 16*data.move(x,y,z,0,2,-1) + 20*data.move(x,y,z,0,1,-1) - 8*data.move(x,y,z,0,0,-1) - 3*data.move(x,y,z,0,3,0) + 12*data.move(x,y,z,0,2,0)
                - 15*data.move(x,y,z,0,1,0) + 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(x,y,z,0,-2,3) + 4*data.move(x,y,z,0,-2,2) - 5*data.move(x,y,z,0,-2,1) + 2*data.move(x,y,z,0,-2,0) + 4*data.move(x,y,z,0,-1,3)
                - 16*data.move(x,y,z,0,-1,2) + 20*data.move(x,y,z,0,-1,1) - 8*data.move(x,y,z,0,-1,0) - 3*data.move(x,y,z,0,0,3) + 12*data.move(x,y,z,0,0,2)
                - 15*data.move(x,y,z,0,0,1) + 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_central_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data.move(x,y,z,-3,1,0) + 4*data.move(x,y,z,-2,1,0) - 5*data.move(x,y,z,-1,1,0) + 2*data.move(x,y,z,0,1,0) + data.move(x,y,z,-3,-1,0)
                - 4*data.move(x,y,z,-2,-1,0) + 5*data.move(x,y,z,-1,-1,0) - 2*data.move(x,y,z,0,-1,0)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data.move(x,y,z,1,-3,0) + 4*data.move(x,y,z,1,-2,0) - 5*data.move(x,y,z,1,-1,0) + 2*data.move(x,y,z,1,0,0) + data.move(x,y,z,-1,-3,0)
                - 4*data.move(x,y,z,-1,-2,0) + 5*data.move(x,y,z,-1,-1,0) - 2*data.move(x,y,z,-1,0,0)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data.move(x,y,z,-3,0,1) + 4*data.move(x,y,z,-2,0,1) - 5*data.move(x,y,z,-1,0,1) + 2*data.move(x,y,z,0,0,1) + data.move(x,y,z,-3,0,-1)
                - 4*data.move(x,y,z,-2,0,-1) + 5*data.move(x,y,z,-1,0,-1) - 2*data.move(x,y,z,0,0,-1)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data.move(x,y,z,1,0,-3) + 4*data.move(x,y,z,1,0,-2) - 5*data.move(x,y,z,1,0,-1) + 2*data.move(x,y,z,1,0,0) + data.move(x,y,z,-1,0,-3)
                - 4*data.move(x,y,z,-1,0,-2) + 5*data.move(x,y,z,-1,0,-1) - 2*data.move(x,y,z,-1,0,0)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data.move(x,y,z,0,-3,1) + 4*data.move(x,y,z,0,-2,1) - 5*data.move(x,y,z,0,-1,1) + 2*data.move(x,y,z,0,0,1) + data.move(x,y,z,0,-3,-1)
                - 4*data.move(x,y,z,0,-2,-1) + 5*data.move(x,y,z,0,-1,-1) - 2*data.move(x,y,z,0,0,-1)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data.move(x,y,z,0,1,-3) + 4*data.move(x,y,z,0,1,-2) - 5*data.move(x,y,z,0,1,-1) + 2*data.move(x,y,z,0,1,0) + data.move(x,y,z,0,-1,-3)
                - 4*data.move(x,y,z,0,-1,-2) + 5*data.move(x,y,z,0,-1,-1) - 2*data.move(x,y,z,0,-1,0)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_forward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz)noexcept  {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data.move(x,y,z,-3,2,0) - 4*data.move(x,y,z,-2,2,0) + 5*data.move(x,y,z,-1,2,0) - 2*data.move(x,y,z,0,2,0) - 4*data.move(x,y,z,-3,1,0)
                + 16*data.move(x,y,z,-2,1,0) - 20*data.move(x,y,z,-1,1,0) + 8*data.move(x,y,z,0,1,0) + 3*data.move(x,y,z,-3,0,0) - 12*data.move(x,y,z,-2,0,0)
                + 15*data.move(x,y,z,-1,0,0) - 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data.move(x,y,z,2,-3,0) - 4*data.move(x,y,z,2,-2,0) + 5*data.move(x,y,z,2,-1,0) - 2*data.move(x,y,z,2,0,0) - 4*data.move(x,y,z,1,-3,0)
                + 16*data.move(x,y,z,1,-2,0) - 20*data.move(x,y,z,1,-1,0) + 8*data.move(x,y,z,1,0,0) + 3*data.move(x,y,z,0,-3,0) - 12*data.move(x,y,z,0,-2,0)
                + 15*data.move(x,y,z,0,-1,0) - 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data.move(x,y,z,-3,0,2) - 4*data.move(x,y,z,-2,0,2) + 5*data.move(x,y,z,-1,0,2) - 2*data.move(x,y,z,0,0,2) - 4*data.move(x,y,z,-3,0,1)
                + 16*data.move(x,y,z,-2,0,1) - 20*data.move(x,y,z,-1,0,1) + 8*data.move(x,y,z,0,0,1) + 3*data.move(x,y,z,-3,0,0) - 12*data.move(x,y,z,-2,0,0)
                + 15*data.move(x,y,z,-1,0,0) - 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data.move(x,y,z,2,0,-3) - 4*data.move(x,y,z,2,0,-2) + 5*data.move(x,y,z,2,0,-1) - 2*data.move(x,y,z,2,0,0) - 4*data.move(x,y,z,1,0,-3)
                + 16*data.move(x,y,z,1,0,-2) - 20*data.move(x,y,z,1,0,-1) + 8*data.move(x,y,z,1,0,0) + 3*data.move(x,y,z,0,0,-3) - 12*data.move(x,y,z,0,0,-2)
                + 15*data.move(x,y,z,0,0,-1) - 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data.move(x,y,z,0,-3,2) - 4*data.move(x,y,z,0,-2,2) + 5*data.move(x,y,z,0,-1,2) - 2*data.move(x,y,z,0,0,2) - 4*data.move(x,y,z,0,-3,1)
                + 16*data.move(x,y,z,0,-2,1) - 20*data.move(x,y,z,0,-1,1) + 8*data.move(x,y,z,0,0,1) + 3*data.move(x,y,z,0,-3,0) - 12*data.move(x,y,z,0,-2,0)
                + 15*data.move(x,y,z,0,-1,0) - 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data.move(x,y,z,0,2,-3) - 4*data.move(x,y,z,0,2,-2) + 5*data.move(x,y,z,0,2,-1) - 2*data.move(x,y,z,0,2,0) - 4*data.move(x,y,z,0,1,-3)
                + 16*data.move(x,y,z,0,1,-2) - 20*data.move(x,y,z,0,1,-1) + 8*data.move(x,y,z,0,1,0) + 3*data.move(x,y,z,0,0,-3) - 12*data.move(x,y,z,0,0,-2)
                + 15*data.move(x,y,z,0,0,-1) - 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}



#endif //CODE_FINITE_DIFFERENCE_HPP
