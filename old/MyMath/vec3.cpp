//
// Created by jacob on 31/7/21.
//

#include "vec3.hpp"
#include <cmath>

 vec3::vec3() : e{0, 0, 0} {}
 vec3::vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}
 vec3::vec3(double E) : e{E, E, E} {}


[[nodiscard]]  double vec3::x() const { return e[0]; }
[[nodiscard]]  double vec3::y() const { return e[1]; }
[[nodiscard]]  double vec3::z() const { return e[2]; }

 vec3 vec3::operator-() const { return vec3(-e[0], -e[1], -e[2]); }
 double vec3::operator[](unsigned i) const { return e[i]; }
 double &vec3::operator[](unsigned i) { return e[i]; }

vec3 &vec3::operator+=(const vec3 &v) {
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];
    return *this;
}

vec3 &vec3::operator*=(const double t) {
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
}

vec3 &vec3::operator/=(const double t) {
    return *this *= 1 / t;
}

 bool vec3::operator==(const vec3 &v) {
    return (e[0] == v.e[0]) && (e[1] == v.e[1]) && (e[2] == v.e[2]);
}

[[nodiscard]]  double vec3::length() const {
    return sqrt(length_squared());
}

[[nodiscard]]  double vec3::length_squared() const {
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
}




//vec3 Utility Functions
std::ostream& operator << (std::ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

 vec3 operator + (const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

 vec3 operator + (const double t, const vec3 &v) {
    return vec3(t+v.e[0], t+v.e[1], t+v.e[2]);
}

 vec3 operator + (const vec3 &v, const double t) {
    return t + v;
}

 vec3 operator - (const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

 vec3 operator * (const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

 vec3 operator * (const double t, const vec3 &v) {
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

 vec3 operator * (const vec3 &v, const double t) {
    return t * v;
}

 vec3 operator / (const vec3 &v, const double t) {
    return (1/t) * v;
}

 double dot(const vec3 &u, const vec3 &v) {
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

 vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3( u.e[1] * v.e[2] - u.e[2] * v.e[1],  u.e[2] * v.e[0] - u.e[0] * v.e[2],  u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

 vec3 unit_vector(const vec3& v) {
    return v / v.length();
}

double quadruple(const vec3& A, const vec3& B) {
    return sqrt(A.length_squared()*B.length_squared() - dot(A,B)*dot(A,B));
}

