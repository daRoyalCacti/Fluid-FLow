//
// Created by jacob on 31/7/21.
//

#ifndef CODE_VEC3_HPP
#define CODE_VEC3_HPP

#include <iostream>

#include <cmath>


struct vec3 {
    double e[3];
    constexpr vec3() : e{0, 0, 0} {}
    constexpr vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}
    explicit constexpr vec3(double E) : e{E, E, E} {}


    [[nodiscard]] constexpr  double x() const { return e[0]; }
    [[nodiscard]] constexpr  double y() const { return e[1]; }
    [[nodiscard]] constexpr  double z() const { return e[2]; }

    constexpr vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    constexpr double operator[](unsigned i) const { return e[i]; }
    constexpr double &operator[](unsigned i) { return e[i]; }

    vec3 &operator+=(const vec3 &v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3 &operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    //for testing
    /*vec3 &operator=(const double d) {
        e[0] = d;
        e[1] = d;
        e[2] = d;
        return *this;
    }*/


    vec3 &operator/=(const double t) {
        return *this *= 1 / t;
    }

    bool operator==(const vec3 &v) {
        return (e[0] == v.e[0]) && (e[1] == v.e[1]) && (e[2] == v.e[2]);
    }

    [[nodiscard]] constexpr  double length() const {
        return sqrt(length_squared());
    }

    [[nodiscard]] constexpr  double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }
};

/*
#include <Eigen/Dense>

//testing
// - https://stackoverflow.com/questions/21112148/specialization-of-member-function-template-after-instantiation-error-and-order
namespace Eigen {

    template<> struct NumTraits<vec3> : GenericNumTraits<vec3>
    {
        typedef vec3 Real;
        typedef vec3 NonInteger;
        typedef vec3 Nested;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 3,
            MulCost = 3
        };
    };

}
 */



//vec3 Utility Functions

std::ostream& operator << (std::ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

constexpr vec3 operator + (const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

constexpr vec3 operator + (const double t, const vec3 &v) {
    return vec3(t+v.e[0], t+v.e[1], t+v.e[2]);
}

constexpr vec3 operator + (const vec3 &v, const double t) {
    return t + v;
}

constexpr vec3 operator - (const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

constexpr vec3 operator * (const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

constexpr vec3 operator * (const double t, const vec3 &v) {
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

constexpr vec3 operator * (const vec3 &v, const double t) {
    return t * v;
}

constexpr vec3 operator / (const vec3 &v, const double t) {
    return (1/t) * v;
}

constexpr double dot(const vec3 &u, const vec3 &v) {
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

constexpr vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3( u.e[1] * v.e[2] - u.e[2] * v.e[1],  u.e[2] * v.e[0] - u.e[0] * v.e[2],  u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

constexpr vec3 unit_vector(const vec3& v) {
    return v / v.length();
}

//computes the quadruple product
//https://mathworld.wolfram.com/VectorQuadrupleProduct.html
constexpr double quadruple(const vec3& A, const vec3& B) {
    return sqrt(A.length_squared()*B.length_squared() - dot(A,B)*dot(A,B));
}


constexpr vec3 tan(const vec3& v) {
    return vec3(tan(v.x()), tan(v.y()), tan(v.z()));
}

constexpr vec3 sin(const vec3& v) {
    return vec3(sin(v.x()), sin(v.y()), sin(v.z()));
}

constexpr vec3 cos(const vec3& v) {
    return vec3(cos(v.x()), cos(v.y()), cos(v.z()));
}

constexpr vec3 sqrt(const vec3& v) {
    return vec3(sqrt(v.x()), sqrt(v.y()), sqrt(v.z()));
}

constexpr vec3 pow(const vec3& v, const double p) {
    return vec3(pow(v.x(), p), pow(v.y(), p), pow(v.z(),p));
}

constexpr vec3 abs2(const vec3& v) {
    return vec3( v.x()*v.x(), v.y()*v.y(), v.z()*v.z());
}

#endif //CODE_VEC3_HPP
