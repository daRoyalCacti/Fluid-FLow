//
// Created by jacob on 24/8/21.
//

#ifndef CODE_TRIANGLE_HPP
#define CODE_TRIANGLE_HPP

#include <utility>
#include "ray.hpp"


struct triangle {
    const vec3 *vertex0, *vertex1, *vertex2;	//position of vertex
    vec3 v0, v1;	//edges of the triangle
    //precomputed quantities to find the uv coordinates
    //https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    double d00{}, d01{}, d11{};	//helpful quantities for finding texture coords

    const vec3 *normal0, *normal1, *normal2;	//vertex normals

    triangle() = delete;


     triangle(const vec3 *vec0, const vec3 *vec1, const vec3 *vec2, const vec3 *n0, const vec3 *n1, const vec3 *n2)
              : vertex0(vec0), vertex1(vec1), vertex2(vec2),  normal0(n0), normal1(n1), normal2(n2) {
         update();
     }


    bool hit_time(const ray& r, double t_min, double t_max, double& hit_time);

     //call whenever the pointers to the data update
     void update() {
         v0 = *vertex1 - *vertex0;
         v1 = *vertex2 - *vertex0;
         d00 = dot(v0, v0)/(dot(v0, v0) * dot(v1, v1) - dot(v0, v1) * dot(v0, v1));
         d01 = dot(v0, v1)/(dot(v0, v0) * dot(v1, v1) - dot(v0, v1) * dot(v0, v1));
         d11 = dot(v1, v1)/(dot(v0, v0) * dot(v1, v1) - dot(v0, v1) * dot(v0, v1));
     }

    inline void barycentric_coords(const vec3 &p, double& Bary0, double& Bary1, double &Bary2) const {
        //find the uv coordinates by interpolating using barycentric coordinates
        //https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
        const vec3 v2 = p - vertex0;
        const double d20 = dot(v2, v0);
        const double d21 = dot(v2, v1);

        //Bary0 = (d11 * d20 - d01 * d21) * invDenom;
        //Bary1 = (d00 * d21 - d01 * d20) * invDenom;
        Bary0 = (d11 * d20 - d01 * d21);    //invDenom now integrated into the definitions of d11, d01 and d00
        Bary1 = (d00 * d21 - d01 * d20);
        Bary2 = 1.0 - Bary0 - Bary1;
    }

    template <typename T>
    static inline T barycentric_interp(const T &interp0, const T &interp1, const T &interp2, const double Bary0, const double Bary1, const double Bary2) {
        //https://computergraphics.stackexchange.com/questions/1866/how-to-map-square-texture-to-triangle
        return Bary2*interp0 + Bary0*interp1 + Bary1*interp2;
    }

    void get_normals(const vec3& point) {
        double Bary0, Bary1, Bary2;
        barycentric_coords(point, Bary0, Bary1, Bary2);
        vec3 temp_norm_res = barycentric_interp(normal0, normal1, normal2, Bary0, Bary1, Bary2);
    }

};

bool triangle::hit_time(const ray& r, const double t_min, const double t_max, double& hit_time) {
    //using the Moller-Trumbore intersection algorithm
    //https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    constexpr double epsilon = 0.0000001;
    vec3 h, s, q;
    double a, f, u, v;


    h = cross(r.dir, v1);
    a = dot(v0, h);

    if (a > -epsilon && a < epsilon)	//ray is parallel to triangle
        return false;

    f = 1.0f / a;
    s = r.orig - vertex0;
    u = f * dot(s, h);

    if (u < 0.0f || u > 1.0f)
        return false;

    q = cross(s, v0);
    v = f * dot(r.dir, q);

    if (v < 0.0f | u + v > 1.0f)
        return false;

    //computing the time of intersection
    const double t = f * dot(v1, q);

    if (t < t_min || t > t_max || t < epsilon)	//time of intersection falls outside of time range considered
        return false;

    hit_time = t;

    return true;
}


#endif //CODE_TRIANGLE_HPP
