//
// Created by jacob on 31/7/21.
//

#include "ops.hpp"

vec3 operator * (const Eigen::Matrix4d &M, const vec3 &v) {
    const Eigen::Vector4d v_E = Eigen::Vector4d(v.x(), v.y(), v.z(), 1);
    const Eigen::Vector4d R = M*v_E;
    return vec3(R.x(), R.y(), R.z());
}

//http://paulbourke.net/geometry/rotate/
vec3 rotate(const vec3&  orig, const vec3& dir, const vec3& point, const double theta) {
    const auto x1 = orig.x();
    const auto y1 = orig.y();
    const auto z1 = orig.z();

    const auto A = dir.x();
    const auto B = dir.y();
    const auto C = dir.z();

    /*const auto x2 = x1 + A;
    const auto y2 = y1 + B;
    const auto z2 = z1 + C;*/

    const auto L = dir.length();
    const auto V = sqrt(B*B + C*C);

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

    const bool V_not_zero = (V > 0.000001);
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
    Ry << V/L, 0, -A/L, 0,
           0,   1,  0,   0,
           A/L, 0,  V/L, 0,
           0,   0,  0,   1;

    Eigen::Matrix4d IRy;
    IRy << V/L, 0, A/L, 0,
           0,   1, 0,   0,
          -A/L, 0, V/L, 0,
           0,   0, 0,   1;

    Eigen::Matrix4d Rz;
    Rz << cos(theta), -sin(theta), 0, 0,
          sin(theta),  cos(theta), 0, 0,
          0,           0,          1, 0,
          0,           0,          0, 1;

    const Eigen::Matrix4d Q = IT * IRx * IRy * Rz * Ry * Rx * T;

    /*std::cout << "================================================\n";
    std::cout << V << "\t" << C << "\t" << V_not_zero << "\n";
    std::cout << IT << "\n\n";
    std::cout << IRx << "\n\n";
    std::cout << IRy << "\n\n";
    std::cout << Rz << "\n\n";
    std::cout << Ry << "\n\n";
    std::cout << Rx << "\n\n";
    std::cout << T << "\n\n";
    std::cout << Q << "\n\n";
    std::cout << "================================================\n";*/

    return Q * point;
}

c_line::c_line() = default;
c_line::c_line(const vec3& x1_, const vec3& x2_) : x1(x1_), x2(x2_), x2_l_x1(x2-x1), dist_denom_inv(1/(x2_-x1_).length()) {}

//https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
double c_line::distance(const vec3& x0) const {
    const auto num = quadruple(x2_l_x1, x1-x0);
    return abs(num) * dist_denom_inv;
}

