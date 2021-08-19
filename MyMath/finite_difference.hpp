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
inline constexpr auto central_difference_1st(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data(x+1,y,z) - data(x-1,y,z) ) / (2*dx);
    } else if constexpr(axis == 1) {
        return ( data(x,y+1,z) - data(x,y-1,z) ) / (2*dy);
    } else if constexpr(axis == 2) {
        return ( data(x,y,z+1) - data(x,y,z-1) ) / (2*dz);
    }
}

template <unsigned axis, typename T>
inline constexpr auto forward_difference_1st(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( -data(x+2,y,z) + 4*data(x+1,y,z) - 3*data(x,y,z)) / (2*dx);
    } else if constexpr(axis == 1) {
        return ( -data(x,y+2,z) + 4*data(x,y+1,z) - 3*data(x,y,z)) / (2*dy);
    } else if constexpr(axis == 2) {
        return ( -data(x,y,z+2) + 4*data(x,y,z+1) - 3*data(x,y,z)) / (2*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto backward_difference_1st(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( 3*data(x,y,z) - 4*data(x-1,y,z) + data(x-2,y,z) ) / (2*dx);
    } else if constexpr(axis == 1) {
        return ( 3*data(x,y,z) - 4*data(x,y-1,z) + data(x,y-2,z) ) / (2*dy);
    } else if constexpr(axis == 2) {
        return ( 3*data(x,y,z) - 4*data(x,y,z-1) + data(x,y,z-2) ) / (2*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto central_difference_2nd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data(x+1,y,z) - 2* data(x,y,z) + data(x-1,y,z) ) / (dx*dx);
    } else if constexpr(axis == 1) {
        return  ( data(x,y+1,z) - 2* data(x,y,z) + data(x,y-1,z) ) / (dy*dy);
    } else if constexpr(axis == 2) {
        return  ( data(x,y,z+1) - 2* data(x,y,z) + data(x,y,z-1) ) / (dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto forward_difference_2nd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (-data(x+3,y,z)+4*data(x+2,y,z) - 5*data(x+1,y,z)+2*data(x,y,z)) / (dx*dx);
    } else if constexpr(axis == 1) {
        return (-data(x,y+3,z)+4*data(x,y+2,z) - 5*data(x,y+1,z)+2*data(x,y,z)) / (dy*dy);
    } else if constexpr(axis == 2) {
        return (-data(x,y,z+3)+4*data(x,y,z+2) - 5*data(x,y,z+1)+2*data(x,y,z)) / (dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto backward_difference_2nd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( 2*data(x,y,z) -5*data(x-1,y,z)+4*data(x-2,y,z) - data(x-3,y,z) ) / (dx*dx);
    } else if constexpr(axis == 1) {
        return ( 2*data(x,y,z) -5*data(x,y-1,z)+4*data(x,y-2,z) - data(x,y-3,z) ) / (dy*dy);
    } else if constexpr(axis == 2) {
        return ( 2*data(x,y,z) -5*data(x,y,z-1)+4*data(x,y,z-2) - data(x,y,z-3) ) / (dz*dz);
    }
}

template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return central_difference_2nd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data(x+1,y+1,z) - data(x-1,y+1,z) - data(x+1,y-1,z) + data(x-1,y-1,z))  / (4*dx*dy);
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data(x+1,y,z+1) - data(x-1,y,z+1) - data(x+1,y,z-1) + data(x-1,y,z-1))  / (4*dx*dz);
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data(x,y+1,z+1) - data(x,y-1,z+1) - data(x,y+1,z-1) + data(x,y-1,z-1))  / (4*dy*dz);
    }
}

template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return forward_difference_2nd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data(x+2,y+2,z) - 4*data(x+1,y+2,z) + 3*data(x,y+2,z) - 4*data(x+2,y+1,z)+16*data(x+1,y+1,z)
                    -12*data(x,y+1,z) + 3*data(x+2,y,z) - 12*data(x+1,y,z) + 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data(x+2,y,z+2) - 4*data(x+1,y,z+2) + 3*data(x,y,z+2) - 4*data(x+2,y,z+1)+16*data(x+1,y,z+1)
                 -12*data(x,y,z+1) + 3*data(x+2,y,z) - 12*data(x+1,y,z) + 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data(x,y+2,z+2) - 4*data(x,y+1,z+2) + 3*data(x,y,z+2) - 4*data(x,y+2,z+1)+16*data(x,y+1,z+1)
                 -12*data(x,y,z+1) + 3*data(x,y+2,z) - 12*data(x,y+1,z) + 9*data(x,y,z)) / (4*dy*dz);
    }
}


template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return backward_difference_2nd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) ) {  //d^/dxdy
        return ( data(x-2,y-2,z) - 4*data(x-1,y-2,z) + 3*data(x,y-2,z) - 4*data(x-2,y-1,z)+16*data(x-1,y-1,z)
                 -12*data(x,y-1,z) + 3*data(x-2,y,z) - 12*data(x-1,y,z) + 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr( (axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0) ) {  //d^2/dxdz
        return ( data(x-2,y,z-2) - 4*data(x-1,y,z-2) + 3*data(x,y,z-2) - 4*data(x-2,y,z-1)+16*data(x-1,y,z-1)
                 -12*data(x,y,z-1) + 3*data(x-2,y,z) - 12*data(x-1,y,z) + 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr( (axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) { //d^2/dydz
        return ( data(x,y-2,z-2) - 4*data(x,y-1,z-2) + 3*data(x,y,z-2) - 4*data(x,y-2,z-1)+16*data(x,y-1,z-1)
                 -12*data(x,y,z-1) + 3*data(x,y-2,z) - 12*data(x,y-1,z) + 9*data(x,y,z)) / (4*dy*dz);
    }
}



template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_forward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return (-data(x+2,y+1,z) + 4*data(x+1,y+1,z) - 3*data(x,y+1,z) + data(x+2,y-1,z) - 4*data(x+1,y-1,z)+3*data(x,y-1,z)) / (4*dx*dy);
    } else if constexpr(axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return (-data(x+1,y+2,z) + 4*data(x+1,y+1,z) - 3*data(x+1,y,z) + data(x-1,y+2,z) - 4*data(x-1,y+1,z)+3*data(x-1,y,z)) / (4*dx*dy);
    } else if constexpr( axis1 == 2 && axis2 == 0) {//d^2/dxdz
        return (-data(x+2,y,z+1) + 4*data(x+1,y,z+1) - 3*data(x,y,z+1) + data(x+2,y,z-1) - 4*data(x+1,y,z-1)+3*data(x,y,z-1)) / (4*dx*dz);
    } else if constexpr(axis1 == 0 && axis2 == 2 ) {  //d^2/dxdz
        return (-data(x+1,y,z+2) + 4*data(x+1,y,z+1) - 3*data(x+1,y,z) + data(x-1,y,z+2) - 4*data(x-1,y,z+1)+3*data(x-1,y,z)) / (4*dx*dz);
    } else if constexpr( axis1 == 2 && axis2 == 1) {//d^2/dydz
        return (-data(x,y+2,z+1) + 4*data(x,y+1,z+1) - 3*data(x,y,z+1) + data(x,y+2,z-1) - 4*data(x,y+1,z-1)+3*data(x,y,z-1)) / (4*dy*dz);
    } else if constexpr(axis1 == 1 && axis2 == 2) { //d^2/dydz
        return (-data(x,y+1,z+2) + 4*data(x,y+1,z+1) - 3*data(x,y+1,z) + data(x,y-1,z+2) - 4*data(x,y-1,z+1)+3*data(x,y-1,z)) / (4*dy*dz);
    }
}


template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_backward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return (data(x-2,y+1,z) - 4*data(x-1,y+1,z) + 3*data(x,y+1,z) - data(x-2,y-1,z) + 4*data(x-1,y-1,z)-3*data(x,y-1,z)) / (4*dx*dy);
    } else if constexpr(axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return (data(x+1,y-2,z) - 4*data(x+1,y-1,z) + 3*data(x+1,y,z) - data(x-1,y-2,z) + 4*data(x-1,y-1,z)-3*data(x-1,y,z)) / (4*dx*dy);
    } else if constexpr( axis1 == 2 && axis2 == 0) {//d^2/dxdz
        return (data(x-2,y,z+1) - 4*data(x-1,y,z+1) + 3*data(x,y,z+1) - data(x-2,y,z-1) + 4*data(x-1,y,z-1)-3*data(x,y,z-1)) / (4*dx*dz);
    } else if constexpr(axis1 == 0 && axis2 == 2 ) {  //d^2/dxdz
        return (data(x+1,y,z-2) - 4*data(x+1,y,z-1) + 3*data(x+1,y,z) - data(x-1,y,z-2) + 4*data(x-1,y,z-1)-3*data(x-1,y,z)) / (4*dx*dz);
    } else if constexpr( axis1 == 2 && axis2 == 1) {//d^2/dydz
        return (data(x,y-2,z+1) - 4*data(x,y-1,z+1) + 3*data(x,y,z+1) - data(x,y-2,z-1) + 4*data(x,y-1,z-1)-3*data(x,y,z-1)) / (4*dy*dz);
    } else if constexpr(axis1 == 1 && axis2 == 2) { //d^2/dydz
        return (data(x,y+1,z-2) - 4*data(x,y+1,z-1) + 3*data(x,y+1,z) - data(x,y-1,z-2) + 4*data(x,y-1,z-1)-3*data(x,y-1,z)) / (4*dy*dz);
    }
}


template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_backward_difference_2nd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^2/dxdy
        return ( -data(x+2,y-2,z) + 4*data(x+1,y-2,z) - 3*data(x,y-2,z) + 4*data(x+2,y-1,z)-16*data(x+1,y-1,z)
                 +12*data(x,y-1,z) - 3*data(x+2,y,z) + 12*data(x+1,y,z) - 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^2/dxdy
        return ( -data(x-2,y+2,z) + 4*data(x-1,y+2,z) - 3*data(x,y+2,z) + 4*data(x-2,y+1,z)-16*data(x-1,y+1,z)
                 +12*data(x,y+1,z) - 3*data(x-2,y,z) + 12*data(x-1,y,z) - 9*data(x,y,z)) / (4*dx*dy);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^2/dxdz
        return ( -data(x+2,y,z-2) + 4*data(x+1,y,z-2) - 3*data(x,y,z-2) + 4*data(x+2,y,z-1)-16*data(x+1,y,z-1)
                 +12*data(x,y,z-1) - 3*data(x+2,y,z) + 12*data(x+1,y,z) - 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^2/dxdz
        return ( -data(x-2,y,z+2) + 4*data(x-1,y,z+2) - 3*data(x,y,z+2) + 4*data(x-2,y,z+1)-16*data(x-1,y,z+1)
                 +12*data(x,y,z+1) - 3*data(x-2,y,z) + 12*data(x-1,y,z) - 9*data(x,y,z)) / (4*dx*dz);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^2/dydz
        return ( -data(x,y+2,z-2) + 4*data(x,y+1,z-2) - 3*data(x,y,z-2) + 4*data(x,y+2,z-1)-16*data(x,y+1,z-1)
                 +12*data(x,y,z-1) - 3*data(x,y+2,z) + 12*data(x,y+1,z) - 9*data(x,y,z)) / (4*dy*dz);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^2/dydz
        return ( -data(x,y-2,z+2) + 4*data(x,y-1,z+2) - 3*data(x,y,z+2) + 4*data(x,y-2,z+1)-16*data(x,y-1,z+1)
                 +12*data(x,y,z+1) - 3*data(x,y-2,z) + 12*data(x,y-1,z) - 9*data(x,y,z)) / (4*dy*dz);
    }
}



template <unsigned axis, typename T>
inline constexpr auto central_difference_3rd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return ( data(x+2,y,z) - 2*data(x+1,y,z) + 2*data(x-1,y,z) - data(x-2,y,z)) / (2*dx*dx*dx);
    } else if constexpr(axis == 1) {
        return ( data(x,y+2,z) - 2*data(x,y+1,z) + 2*data(x,y-1,z) - data(x,y-2,z)) / (2*dy*dy*dy);
    } else if constexpr(axis == 2) {
        return ( data(x,y,z+2) - 2*data(x,y,z+1) + 2*data(x,y,z-1) - data(x,y,z-2)) / (2*dz*dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto forward_difference_3rd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (-3*data(x+4,y,z)+14*data(x+3,y,z) - 24*data(x+2,y,z)+ 18*data(x+1,y,z) - 5*data(x,y,z)) / (2*dx*dx*dx);
    } else if constexpr(axis == 1) {
        return (-3*data(x,y+4,z)+14*data(x,y+3,z) - 24*data(x,y+2,z)+ 18*data(x,y+1,z) - 5*data(x,y,z)) / (2*dy*dy*dy);
    } else if constexpr(axis == 2) {
        return (-3*data(x,y,z+4)+14*data(x,y,z+3) - 24*data(x,y,z+2)+ 18*data(x,y,z+1) - 5*data(x,y,z)) / (2*dz*dz*dz);
    }
}


template <unsigned axis, typename T>
inline constexpr auto backward_difference_3rd(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis < 3, "Axis should be 0,1,2");
    if constexpr (axis == 0) {
        return (5*data(x,y,z) - 18*data(x-1,y,z) + 24*data(x-2,y,z) - 14*data(x-3,y,z) + 3*data(x-4,y,z)) / (2*dx*dx*dx);
    } else if constexpr(axis == 1) {
        return (5*data(x,y,z) - 18*data(x,y-1,z) + 24*data(x,y-2,z) - 14*data(x,y-3,z) + 3*data(x,y-4,z)) / (2*dy*dy*dy);
    } else if constexpr(axis == 2) {
        return (5*data(x,y,z) - 18*data(x,y,z-1) + 24*data(x,y,z-2) - 14*data(x,y,z-3) + 3*data(x,y,z-4)) / (2*dz*dz*dz);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return central_difference_3rd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data(x+1,y+1,z) - 2*data(x,y+1,z) + data(x-1,y+1,z) - data(x+1,y-1,z)+2*data(x,y-1,z)-data(x-1,y-1,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data(x+1,y+1,z) - 2*data(x+1,y,z) + data(x+1,y-1,z) - data(x-1,y+1,z)+2*data(x-1,y,z)-data(x-1,y-1,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data(x+1,y,z+1) - 2*data(x,y,z+1) + data(x-1,y,z+1) - data(x+1,y,z-1)+2*data(x,y,z-1)-data(x-1,y,z-1)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data(x+1,y,z+1) - 2*data(x+1,y,z) + data(x+1,y,z-1) - data(x-1,y,z+1)+2*data(x-1,y,z)-data(x-1,y,z-1)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data(x,y+1,z+1) - 2*data(x,y,z+1) + data(x,y-1,z+1) - data(x,y+1,z-1)+2*data(x,y,z-1)-data(x,y-1,z-1)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data(x,y+1,z+1) - 2*data(x,y+1,z) + data(x,y+1,z-1) - data(x,y-1,z+1)+2*data(x,y-1,z)-data(x,y-1,z-1)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return forward_difference_3rd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data(x+3,y+2,z) - 4*data(x+2,y+2,z) + 5*data(x+1,y+2,z) - 2*data(x,y+2,z) - 4*data(x+3,y+1,z)
                 + 16*data(x+2,y+1,z) - 20*data(x+1,y+1,z) + 8*data(x,y+1,z) + 3*data(x+3,y,z) - 12*data(x+2,y,z)
                 + 15*data(x+1,y,z) - 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data(x+2,y+3,z) - 4*data(x+2,y+2,z) + 5*data(x+2,y+1,z) - 2*data(x+2,y,z) - 4*data(x+1,y+3,z)
                + 16*data(x+1,y+2,z) - 20*data(x+1,y+1,z) + 8*data(x+1,y,z) + 3*data(x,y+3,z) - 12*data(x,y+2,z)
                + 15*data(x,y+1,z) - 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data(x+3,y,z+2) - 4*data(x+2,y,z+2) + 5*data(x+1,y,z+2) - 2*data(x,y,z+2) - 4*data(x+3,y,z+1)
                + 16*data(x+2,y,z+1) - 20*data(x+1,y,z+1) + 8*data(x,y,z+1) + 3*data(x+3,y,z) - 12*data(x+2,y,z)
                + 15*data(x+1,y,z) - 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data(x+2,y,z+3) - 4*data(x+2,y,z+2) + 5*data(x+2,y,z+1) - 2*data(x+2,y,z) - 4*data(x+1,y,z+3)
                + 16*data(x+1,y,z+2) - 20*data(x+1,y,z+1) + 8*data(x+1,y,z) + 3*data(x,y,z+3) - 12*data(x,y,z+2)
                + 15*data(x,y,z+1) - 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data(x,y+3,z+2) - 4*data(x,y+2,z+2) + 5*data(x,y+1,z+2) - 2*data(x,y,z+2) - 4*data(x,y+3,z+1)
                + 16*data(x,y+2,z+1) - 20*data(x,y+1,z+1) + 8*data(x,y,z+1) + 3*data(x,y+3,z) - 12*data(x,y+2,z)
                + 15*data(x,y+1,z) - 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data(x,y+2,z+3) - 4*data(x,y+2,z+2) + 5*data(x,y+2,z+1) - 2*data(x,y+2,z) - 4*data(x,y+1,z+3)
                + 16*data(x,y+1,z+2) - 20*data(x,y+1,z+1) + 8*data(x,y+1,z) + 3*data(x,y,z+3) - 12*data(x,y,z+2)
                + 15*data(x,y,z+1) - 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");

    if constexpr(axis1 == axis2) {
        return backward_difference_3rd<axis1>(data, x, y, z, dx, dy, dz);
    } else if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data(x-3,y-2,z) + 4*data(x-2,y-2,z) - 5*data(x-1,y-2,z) + 2*data(x,y-2,z) + 4*data(x-3,y-1,z)
                - 16*data(x-2,y-1,z) + 20*data(x-1,y-1,z) - 8*data(x,y-1,z) - 3*data(x-3,y,z) + 12*data(x-2,y,z)
                - 15*data(x-1,y,z) + 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data(x-2,y-3,z) + 4*data(x-2,y-2,z) - 5*data(x-2,y-1,z) + 2*data(x-2,y,z) + 4*data(x-1,y-3,z)
                - 16*data(x-1,y-2,z) + 20*data(x-1,y-1,z) - 8*data(x-1,y,z) - 3*data(x,y-3,z) + 12*data(x,y-2,z)
                - 15*data(x,y-1,z) + 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data(x-3,y,z-2) + 4*data(x-2,y,z-2) - 5*data(x-1,y,z-2) + 2*data(x,y,z-2) + 4*data(x-3,y,z-1)
                - 16*data(x-2,y,z-1) + 20*data(x-1,y,z-1) - 8*data(x,y,z-1) - 3*data(x-3,y,z) + 12*data(x-2,y,z)
                - 15*data(x-1,y,z) + 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data(x-2,y,z-3) + 4*data(x-2,y,z-2) - 5*data(x-2,y,z-1) + 2*data(x-2,y,z) + 4*data(x-1,y,z-3)
                - 16*data(x-1,y,z-2) + 20*data(x-1,y,z-1) - 8*data(x-1,y,z) - 3*data(x,y,z-3) + 12*data(x,y,z-2)
                - 15*data(x,y,z-1) + 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data(x,y-3,z-2) + 4*data(x,y-2,z-2) - 5*data(x,y-1,z-2) + 2*data(x,y,z-2) + 4*data(x,y-3,z-1)
                - 16*data(x,y-2,z-1) + 20*data(x,y-1,z-1) - 8*data(x,y,z-1) - 3*data(x,y-3,z) + 12*data(x,y-2,z)
                - 15*data(x,y-1,z) + 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data(x,y-2,z-3) + 4*data(x,y-2,z-2) - 5*data(x,y-2,z-1) + 2*data(x,y-2,z) + 4*data(x,y-1,z-3)
                - 16*data(x,y-1,z-2) + 20*data(x,y-1,z-1) - 8*data(x,y-1,z) - 3*data(x,y,z-3) + 12*data(x,y,z-2)
                - 15*data(x,y,z-1) + 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_forward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data(x+1,y+2,z) + 2*data(x,y+2,z) - data(x-1,y+2,z) + 4*data(x+1,y+1,z) - 8*data(x,y+1,z)
                 + 4*data(x-1,y+1,z) - 3*data(x+1,y,z) + 6*data(x,y,z) - 3*data(x-1,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data(x+2,y+1,z) + 2*data(x+2,y,z) - data(x+2,y-1,z) + 4*data(x+1,y+1,z) - 8*data(x+1,y,z)
                + 4*data(x+1,y-1,z) - 3*data(x,y+1,z) + 6*data(x,y,z) - 3*data(x,y-1,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data(x+1,y,z+2) + 2*data(x,y,z+2) - data(x-1,y,z+2) + 4*data(x+1,y,z+1) - 8*data(x,y,z+1)
                + 4*data(x-1,y,z+1) - 3*data(x+1,y,z) + 6*data(x,y,z) - 3*data(x-1,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data(x+2,y,z+1) + 2*data(x+2,y,z) - data(x+2,y,z-1) + 4*data(x+1,y,z+1) - 8*data(x+1,y,z)
                + 4*data(x+1,y,z-1) - 3*data(x,y,z+1) + 6*data(x,y,z) - 3*data(x,y,z-1)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data(x,y+1,z+2) + 2*data(x,y,z+2) - data(x,y-1,z+2) + 4*data(x,y+1,z+1) - 8*data(x,y,z+1)
                + 4*data(x,y-1,z+1) - 3*data(x,y+1,z) + 6*data(x,y,z) - 3*data(x,y-1,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data(x,y+2,z+1) + 2*data(x,y+2,z) - data(x,y+2,z-1) + 4*data(x,y+1,z+1) - 8*data(x,y+1,z)
                + 4*data(x,y+1,z-1) - 3*data(x,y,z+1) + 6*data(x,y,z) - 3*data(x,y,z-1)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto central_backward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data(x+1,y-2,z) - 2*data(x,y-2,z) + data(x-1,y-2,z) - 4*data(x+1,y-1,z) + 8*data(x,y-1,z)
                - 4*data(x-1,y-1,z) + 3*data(x+1,y,z) - 6*data(x,y,z) + 3*data(x-1,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data(x-2,y+1,z) - 2*data(x-2,y,z) + data(x-2,y-1,z) - 4*data(x-1,y+1,z) + 8*data(x-1,y,z)
                - 4*data(x-1,y-1,z) + 3*data(x,y+1,z) - 6*data(x,y,z) + 3*data(x,y-1,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data(x+1,y,z-2) - 2*data(x,y,z-2) + data(x-1,y,z-2) - 4*data(x+1,y,z-1) + 8*data(x,y,z-1)
                - 4*data(x-1,y,z-1) + 3*data(x+1,y,z) - 6*data(x,y,z) + 3*data(x-1,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data(x-2,y,z+1) - 2*data(x-2,y,z) + data(x-2,y,z-1) - 4*data(x-1,y,z+1) + 8*data(x-1,y,z)
                - 4*data(x-1,y,z-1) + 3*data(x,y,z+1) - 6*data(x,y,z) + 3*data(x,y,z-1)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data(x,y+1,z-2) - 2*data(x,y,z-2) + data(x,y-1,z-2) - 4*data(x,y+1,z-1) + 8*data(x,y,z-1)
                - 4*data(x,y-1,z-1) + 3*data(x,y+1,z) - 6*data(x,y,z) + 3*data(x,y-1,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data(x,y-2,z+1) - 2*data(x,y-2,z) + data(x,y-2,z-1) - 4*data(x,y-1,z+1) + 8*data(x,y-1,z)
                - 4*data(x,y-1,z-1) + 3*data(x,y,z+1) - 6*data(x,y,z) + 3*data(x,y,z-1)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_central_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data(x+3,y+1,z) + 4*data(x+2,y+1,z) - 5*data(x+1,y+1,z) + 2*data(x,y+1,z) + data(x+3,y-1,z)
                 - 4*data(x+2,y-1,z) + 5*data(x+1,y-1,z) - 2*data(x,y-1,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data(x+1,y+3,z) + 4*data(x+1,y+2,z) - 5*data(x+1,y+1,z) + 2*data(x+1,y,z) + data(x-1,y+3,z)
                - 4*data(x-1,y+2,z) + 5*data(x-1,y+1,z) - 2*data(x-1,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data(x+3,y,z+1) + 4*data(x+2,y,z+1) - 5*data(x+1,y,z+1) + 2*data(x,y,z+1) + data(x+3,y,z-1)
                - 4*data(x+2,y,z-1) + 5*data(x+1,y,z-1) - 2*data(x,y,z-1)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data(x+1,y,z+3) + 4*data(x+1,y,z+2) - 5*data(x+1,y,z+1) + 2*data(x+1,y,z) + data(x-1,y,z+3)
                - 4*data(x-1,y,z+2) + 5*data(x-1,y,z+1) - 2*data(x-1,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data(x,y+3,z+1) + 4*data(x,y+2,z+1) - 5*data(x,y+1,z+1) + 2*data(x,y,z+1) + data(x,y+3,z-1)
                - 4*data(x,y+2,z-1) + 5*data(x,y+1,z-1) - 2*data(x,y,z-1)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data(x,y+1,z+3) + 4*data(x,y+1,z+2) - 5*data(x,y+1,z+1) + 2*data(x,y+1,z) + data(x,y-1,z+3)
                - 4*data(x,y-1,z+2) + 5*data(x,y-1,z+1) - 2*data(x,y-1,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto forward_backward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data(x+3,y-2,z) + 4*data(x+2,y-2,z) - 5*data(x+1,y-2,z) + 2*data(x,y-2,z) + 4*data(x+3,y-1,z)
         - 16*data(x+2,y-1,z) + 20*data(x+1,y-1,z) - 8*data(x,y-1,z) - 3*data(x+3,y,z) + 12*data(x+2,y,z)
         - 15*data(x+1,y,z) + 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data(x-2,y+3,z) + 4*data(x-2,y+2,z) - 5*data(x-2,y+1,z) + 2*data(x-2,y,z) + 4*data(x-1,y+3,z)
                - 16*data(x-1,y+2,z) + 20*data(x-1,y+1,z) - 8*data(x-1,y,z) - 3*data(x,y+3,z) + 12*data(x,y+2,z)
                - 15*data(x,y+1,z) + 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data(x+3,y,z-2) + 4*data(x+2,y,z-2) - 5*data(x+1,y,z-2) + 2*data(x,y,z-2) + 4*data(x+3,y,z-1)
                - 16*data(x+2,y,z-1) + 20*data(x+1,y,z-1) - 8*data(x,y,z-1) - 3*data(x+3,y,z) + 12*data(x+2,y,z)
                - 15*data(x+1,y,z) + 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data(x-2,y,z+3) + 4*data(x-2,y,z+2) - 5*data(x-2,y,z+1) + 2*data(x-2,y,z) + 4*data(x-1,y,z+3)
                - 16*data(x-1,y,z+2) + 20*data(x-1,y,z+1) - 8*data(x-1,y,z) - 3*data(x,y,z+3) + 12*data(x,y,z+2)
                - 15*data(x,y,z+1) + 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data(x,y+3,z-2) + 4*data(x,y+2,z-2) - 5*data(x,y+1,z-2) + 2*data(x,y,z-2) + 4*data(x,y+3,z-1)
                - 16*data(x,y+2,z-1) + 20*data(x,y+1,z-1) - 8*data(x,y,z-1) - 3*data(x,y+3,z) + 12*data(x,y+2,z)
                - 15*data(x,y+1,z) + 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data(x,y-2,z+3) + 4*data(x,y-2,z+2) - 5*data(x,y-2,z+1) + 2*data(x,y-2,z) + 4*data(x,y-1,z+3)
                - 16*data(x,y-1,z+2) + 20*data(x,y-1,z+1) - 8*data(x,y-1,z) - 3*data(x,y,z+3) + 12*data(x,y,z+2)
                - 15*data(x,y,z+1) + 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}


//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_central_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (-data(x-3,y+1,z) + 4*data(x-2,y+1,z) - 5*data(x-1,y+1,z) + 2*data(x,y+1,z) + data(x-3,y-1,z)
                - 4*data(x-2,y-1,z) + 5*data(x-1,y-1,z) - 2*data(x,y-1,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (-data(x+1,y-3,z) + 4*data(x+1,y-2,z) - 5*data(x+1,y-1,z) + 2*data(x+1,y,z) + data(x-1,y-3,z)
                - 4*data(x-1,y-2,z) + 5*data(x-1,y-1,z) - 2*data(x-1,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (-data(x-3,y,z+1) + 4*data(x-2,y,z+1) - 5*data(x-1,y,z+1) + 2*data(x,y,z+1) + data(x-3,y,z-1)
                - 4*data(x-2,y,z-1) + 5*data(x-1,y,z-1) - 2*data(x,y,z-1)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (-data(x+1,y,z-3) + 4*data(x+1,y,z-2) - 5*data(x+1,y,z-1) + 2*data(x+1,y,z) + data(x-1,y,z-3)
                - 4*data(x-1,y,z-2) + 5*data(x-1,y,z-1) - 2*data(x-1,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (-data(x,y-3,z+1) + 4*data(x,y-2,z+1) - 5*data(x,y-1,z+1) + 2*data(x,y,z+1) + data(x,y-3,z-1)
                - 4*data(x,y-2,z-1) + 5*data(x,y-1,z-1) - 2*data(x,y,z-1)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (-data(x,y+1,z-3) + 4*data(x,y+1,z-2) - 5*data(x,y+1,z-1) + 2*data(x,y+1,z) + data(x,y-1,z-3)
                - 4*data(x,y-1,z-2) + 5*data(x,y-1,z-1) - 2*data(x,y-1,z)) / (2*dz*dz*dy);
    }
}

//first axis is squared
template <unsigned axis1, unsigned axis2, typename T>
inline constexpr auto backward_forward_difference_3rd_mixed(const T& data, const size_t x, const size_t y, const size_t z, const double dx, const double dy, const double dz) {
    static_assert(axis1 < 3, "Axis should be 0,1,2");
    static_assert(axis2 < 3, "Axis should be 0,1,2");
    static_assert(axis1 != axis2, "Cannot have the axes the same");

    if constexpr( axis1 == 0 && axis2 == 1) {  //d^3/dx^2dy
        return (data(x-3,y+2,z) - 4*data(x-2,y+2,z) + 5*data(x-1,y+2,z) - 2*data(x,y+2,z) - 4*data(x-3,y+1,z)
                + 16*data(x-2,y+1,z) - 20*data(x-1,y+1,z) + 8*data(x,y+1,z) + 3*data(x-3,y,z) - 12*data(x-2,y,z)
                + 15*data(x-1,y,z) - 6*data(x,y,z)) / (2*dy*dx*dx);
    } else if constexpr(axis1 == 1 && axis2 == 0) {  //d^3/dy^2dx
        return (data(x+2,y-3,z) - 4*data(x+2,y-2,z) + 5*data(x+2,y-1,z) - 2*data(x+2,y,z) - 4*data(x+1,y-3,z)
                + 16*data(x+1,y-2,z) - 20*data(x+1,y-1,z) + 8*data(x+1,y,z) + 3*data(x,y-3,z) - 12*data(x,y-2,z)
                + 15*data(x,y-1,z) - 6*data(x,y,z)) / (2*dy*dy*dx);
    } else if constexpr( axis1 == 0 && axis2 == 2) {//d^3/dx^2dz
        return (data(x-3,y,z+2) - 4*data(x-2,y,z+2) + 5*data(x-1,y,z+2) - 2*data(x,y,z+2) - 4*data(x-3,y,z+1)
                + 16*data(x-2,y,z+1) - 20*data(x-1,y,z+1) + 8*data(x,y,z+1) + 3*data(x-3,y,z) - 12*data(x-2,y,z)
                + 15*data(x-1,y,z) - 6*data(x,y,z)) / (2*dz*dx*dx);
    } else if constexpr(axis1 == 2 && axis2 == 0 ) {  //d^3/dz^2dx
        return (data(x+2,y,z-3) - 4*data(x+2,y,z-2) + 5*data(x+2,y,z-1) - 2*data(x+2,y,z) - 4*data(x+1,y,z-3)
                + 16*data(x+1,y,z-2) - 20*data(x+1,y,z-1) + 8*data(x+1,y,z) + 3*data(x,y,z-3) - 12*data(x,y,z-2)
                + 15*data(x,y,z-1) - 6*data(x,y,z)) / (2*dz*dz*dx);
    } else if constexpr( axis1 == 1 && axis2 == 2) {//d^3/dy^2dz
        return (data(x,y-3,z+2) - 4*data(x,y-2,z+2) + 5*data(x,y-1,z+2) - 2*data(x,y,z+2) - 4*data(x,y-3,z+1)
                + 16*data(x,y-2,z+1) - 20*data(x,y-1,z+1) + 8*data(x,y,z+1) + 3*data(x,y-3,z) - 12*data(x,y-2,z)
                + 15*data(x,y-1,z) - 6*data(x,y,z)) / (2*dz*dy*dy);
    } else if constexpr(axis1 == 2 && axis2 == 1) { //d^3/dz^2dy
        return (data(x,y+2,z-3) - 4*data(x,y+2,z-2) + 5*data(x,y+2,z-1) - 2*data(x,y+2,z) - 4*data(x,y+1,z-3)
                + 16*data(x,y+1,z-2) - 20*data(x,y+1,z-1) + 8*data(x,y+1,z) + 3*data(x,y,z-3) - 12*data(x,y,z-2)
                + 15*data(x,y,z-1) - 6*data(x,y,z)) / (2*dz*dz*dy);
    }
}



#endif //CODE_FINITE_DIFFERENCE_HPP
