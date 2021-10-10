//
// Created by jacob on 10/10/21.
//

#ifndef CODE_DIST_TO_PLANE_HPP
#define CODE_DIST_TO_PLANE_HPP

#include "vec3.hpp"

constexpr double dist_to_plane(const vec3& x, const vec3 &p1, const vec3 &p2, const vec3 &p3) {
    //https://tutorial.math.lamar.edu/classes/calcIII/EqnsOfPlanes.aspx
    const auto v1 = p1 - p2;
    const auto v2 = p3 - p2;
    const auto n = cross(v1, v2);

    const auto a = n.x();
    const auto b = n.y();
    const auto c = n.z();

    const auto d = a*p1.x() + b*p1.y() + c*p1.z();

    //https://mathworld.wolfram.com/Point-PlaneDistance.html
    const auto num = std::abs( a*x.x() + b*x.y() + c*x.z() - d);
    const auto denom = sqrt(a*a + b*b + c*c);

    return num/denom;
}

#endif //CODE_DIST_TO_PLANE_HPP
