//
// Created by jacob on 31/7/21.
//

#ifndef CODE_OPS_HPP
#define CODE_OPS_HPP

#include "vec3.hpp"
#include <Eigen/Dense>

//const line struct
// - used to precompute some values
// - and for computing the distance a point is from a line
struct c_line final {
    const vec3 x1;
    const vec3 x2;
    const vec3 x2_l_x1;
    //variable that stores the denominator of the fraction when computing the distance from the line to a point
    const double dist_denom_inv = 0;    //distance function will return NaN's if this isn't set

    c_line() = default;
    constexpr c_line(const vec3& x1_, const vec3& x2_) noexcept : x1(x1_), x2(x2_), x2_l_x1(x2-x1), dist_denom_inv(1/(x2_-x1_).length()) {
#ifndef NDEBUG
        if (!std::isfinite(dist_denom_inv)) {
            std::cerr << "Cline gave infinite inverse distance\n";
        }
#endif
    }

    //the distance a point x0 is from the line
    //https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    [[nodiscard]] constexpr double distance(const vec3& x0) const noexcept {
        const auto num = quadruple(x2_l_x1, x1-x0);
#ifndef NDEBUG
        const auto ret = abs(num) * dist_denom_inv;
        if (!std::isfinite(ret)) {
            std::cerr << "cline distance trying to return infinite value\n";
            std::cerr << "\tnum : " << num << "\tret : " << ret << "\n";
        }
        return ret;
#else
        return abs(num) * dist_denom_inv;
#endif
    }
};


//cannot multiply by vec3 so this instead multiplies by the vector (v, 1)
vec3 operator * (const Eigen::Matrix4d &M, const vec3 &v) noexcept {
    const Eigen::Vector4d v_E = Eigen::Vector4d(v.x(), v.y(), v.z(), 1);
    const Eigen::Vector4d R = M*v_E;
    return vec3(R.x(), R.y(), R.z());
}

//http://paulbourke.net/geometry/rotate/
//rotates point about an axis defined by orig and dir (using the right hand rule)
vec3 rotate(const vec3&  orig, const vec3& dir, const vec3& point, const double theta) noexcept {
    const auto x1 = orig.x();
    const auto y1 = orig.y();
    const auto z1 = orig.z();

    const auto A = dir.x();
    const auto B = dir.y();
    const auto C = dir.z();

    const auto L = dir.length();
    const auto V = sqrt(B*B + C*C);

    const bool V_not_zero = (V > 0.000001);
    const bool L_not_zero = (L > 0.000001);

    Eigen::Matrix4d T;
    T << 1, 0, 0, -x1,
            0, 1, 0, -y1,
            0, 0, 1, -z1,
            0, 0, 0, 1;

    Eigen::Matrix4d IT;
    IT << 1, 0, 0, x1,
            0, 1, 0, y1,
            0, 0, 1, z1,
            0, 0, 0, 1;


    Eigen::Matrix4d Rx;
    if (V_not_zero) {  //can get division by 0
        Rx << 1, 0, 0, 0,
                0, C / V, -B / V, 0,
                0, B / V, C / V, 0,
                0, 0, 0, 1;
    } else {
        Rx << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
    }

    Eigen::Matrix4d IRx;
    if (V_not_zero) {   //can get division by 0
        IRx << 1, 0, 0, 0,
                0, C / V, B / V, 0,
                0, -B / V, C / V, 0,
                0, 0, 0, 1;
    } else {
        IRx << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
    }

    Eigen::Matrix4d Ry;
    if (L_not_zero) {
        Ry << V/L, 0, -A/L, 0,
                0,   1,  0,   0,
                A/L, 0,  V/L, 0,
                0,   0,  0,   1;
    } else {
        Ry << 1, 0, 0, 0,
                0,   1,  0,   0,
                0, 0,  1, 0,
                0,   0,  0,   1;
    }

    Eigen::Matrix4d IRy;
    if (L_not_zero) {
        IRy << V/L, 0, A/L, 0,
                0,   1, 0,   0,
                -A/L, 0, V/L, 0,
                0,   0, 0,   1;
    } else {
        IRy << 1, 0, 0, 0,
                0,   1, 0,   0,
                0, 0, 1, 0,
                0,   0, 0,   1;
    }

    Eigen::Matrix4d Rz;
    Rz << cos(theta), -sin(theta), 0, 0,
            sin(theta),  cos(theta), 0, 0,
            0,           0,          1, 0,
            0,           0,          0, 1;

    const Eigen::Matrix4d Q = IT * IRx * IRy * Rz * Ry * Rx * T;

    return Q * point;
}






#endif //CODE_OPS_HPP
