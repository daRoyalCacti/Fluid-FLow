//
// Created by jacob on 31/7/21.
//

#ifndef CODE_OPS_HPP
#define CODE_OPS_HPP

#include "vec3.hpp"
#include <Eigen/Dense>

//cannot multiply by vec3 so this instead multiplies by the vector (v, 1)
vec3 operator * (const Eigen::Matrix4d &M, const vec3 &v);   //not the reverse is not defined

//rotates point about an axis defined by orig and dir (using the right hand rule)
vec3 rotate(const vec3&  orig, const vec3& dir, const vec3& point, double theta);

//const line struct
// - used to precompute some values
// - and for computing the distance a point is from a line
struct c_line {
    const vec3 x1;
    const vec3 x2;
    const vec3 x2_l_x1;
    //variable that stores the denominator of the fraction when computing the distance from the line to a point
    const double dist_denom_inv = 0;    //distance function will return NaN's if this isn't set

    c_line();
    c_line(const vec3& x1_, const vec3& x2_);

    //the distance a point x0 is from the line
    [[nodiscard]] double distance(const vec3& x0) const;
};

#endif //CODE_OPS_HPP
