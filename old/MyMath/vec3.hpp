//
// Created by jacob on 31/7/21.
//

#ifndef CODE_VEC3_HPP
#define CODE_VEC3_HPP

#include <iostream>


struct vec3 {
    double e[3];
     vec3();
     vec3(double e0, double e1, double e2);
    explicit  vec3(double E);


    [[nodiscard]]  double x() const;
    [[nodiscard]]  double y() const;
    [[nodiscard]]  double z() const;

     vec3 operator-() const;
     double operator[](unsigned i) const;
     double &operator[](unsigned i);

    vec3 &operator+=(const vec3 &v);

    vec3 &operator*=(double t);

    vec3 &operator/=(double t);

     bool operator==(const vec3 &v);

    [[nodiscard]]  double length() const;

    [[nodiscard]]  double length_squared() const;
};



//vec3 Utility Functions

std::ostream& operator << (std::ostream &out, const vec3 &v);

 vec3 operator + (const vec3 &u, const vec3 &v);
 vec3 operator + (double t, const vec3 &v);
 vec3 operator + (const vec3 &v, double t);

 vec3 operator - (const vec3 &u, const vec3 &v);

 vec3 operator * (const vec3 &u, const vec3 &v);
 vec3 operator * (double t, const vec3 &v);
 vec3 operator * (const vec3 &v, double t);

 vec3 operator / (const vec3 &v, double t);

 double dot(const vec3 &u, const vec3 &v);

 vec3 cross(const vec3 &u, const vec3 &v);

 vec3 unit_vector(const vec3& v);

 //computes the quadruple product
 //https://mathworld.wolfram.com/VectorQuadrupleProduct.html
 double quadruple(const vec3& A, const vec3& B);

#endif //CODE_VEC3_HPP
