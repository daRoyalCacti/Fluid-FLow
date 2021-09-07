//
// Created by jacob on 31/7/21.
//

#ifndef CODE_VEC3_HPP
#define CODE_VEC3_HPP

#include <iostream>

#include <cmath>


struct vec3 final {
    double e[3];
    constexpr vec3() noexcept : e{0, 0, 0} {}
    constexpr vec3(double e0, double e1, double e2) noexcept : e{e0, e1, e2} {}
    explicit constexpr vec3(double E) noexcept : e{E, E, E} {}


    [[nodiscard]] constexpr  double x() const noexcept { return e[0]; }
    [[nodiscard]] constexpr  double y() const noexcept { return e[1]; }
    [[nodiscard]] constexpr  double z() const noexcept { return e[2]; }

    constexpr vec3 operator-() const noexcept { return vec3(-e[0], -e[1], -e[2]); }
    constexpr double operator[](unsigned i) const noexcept { return e[i]; }
    constexpr double &operator[](unsigned i) noexcept { return e[i]; }

    vec3 &operator+=(const vec3 &v) noexcept {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3 &operator*=(const double t) noexcept {
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


    vec3 &operator/=(const double t) noexcept {
        return *this *= 1 / t;
    }

    bool operator==(const vec3 &v) const noexcept {
        return (e[0] == v.e[0]) && (e[1] == v.e[1]) && (e[2] == v.e[2]);
    }

    [[nodiscard]] constexpr  double length() const noexcept {
        return sqrt(length_squared());
    }

    [[nodiscard]] constexpr  double length_squared() const noexcept {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }
};


//vec3 Utility Functions

std::ostream& operator << (std::ostream &out, const vec3 &v) noexcept {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

constexpr vec3 operator + (const vec3 &u, const vec3 &v) noexcept {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

constexpr vec3 operator + (const double t, const vec3 &v) noexcept {
    return vec3(t+v.e[0], t+v.e[1], t+v.e[2]);
}

constexpr vec3 operator + (const vec3 &v, const double t) noexcept {
    return t + v;
}

constexpr vec3 operator - (const vec3 &u, const vec3 &v) noexcept {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

constexpr vec3 operator * (const vec3 &u, const vec3 &v) noexcept {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

constexpr vec3 operator * (const double t, const vec3 &v) noexcept {
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

constexpr vec3 operator * (const vec3 &v, const double t) noexcept {
    return t * v;
}

constexpr vec3 operator / (const vec3 &v, const double t) noexcept {
    return (1/t) * v;
}

constexpr double dot(const vec3 &u, const vec3 &v) noexcept {
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

constexpr vec3 cross(const vec3 &u, const vec3 &v) noexcept {
    return vec3( u.e[1] * v.e[2] - u.e[2] * v.e[1],  u.e[2] * v.e[0] - u.e[0] * v.e[2],  u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

constexpr vec3 unit_vector(const vec3& v) noexcept {
    return v / v.length();
}

//computes the quadruple product
//https://mathworld.wolfram.com/VectorQuadrupleProduct.html
constexpr double quadruple(const vec3& A, const vec3& B) noexcept {
#ifndef NDEBUG
    const auto a =A.length_squared()*B.length_squared() - dot(A,B)*dot(A,B);
    if (a < 0) {
        std::cerr << "quadruple trying to return complex number\n";
    }
    return sqrt(a);
#else
    return sqrt(A.length_squared()*B.length_squared() - dot(A,B)*dot(A,B));
#endif
}


constexpr vec3 tan(const vec3& v) noexcept {
    return vec3(tan(v.x()), tan(v.y()), tan(v.z()));
}

constexpr vec3 sin(const vec3& v) noexcept {
    return vec3(sin(v.x()), sin(v.y()), sin(v.z()));
}

constexpr vec3 cos(const vec3& v) noexcept {
    return vec3(cos(v.x()), cos(v.y()), cos(v.z()));
}

constexpr vec3 sqrt(const vec3& v) noexcept {
    return vec3(sqrt(v.x()), sqrt(v.y()), sqrt(v.z()));
}

constexpr vec3 pow(const vec3& v, const double p) noexcept {
    return vec3(pow(v.x(), p), pow(v.y(), p), pow(v.z(),p));
}

constexpr vec3 abs2(const vec3& v)noexcept  {
    return vec3( v.x()*v.x(), v.y()*v.y(), v.z()*v.z());
}

#endif //CODE_VEC3_HPP
