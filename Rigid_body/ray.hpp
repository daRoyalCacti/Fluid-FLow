//
// Created by jacob on 24/8/21.
//

#ifndef CODE_RAY_HPP
#define CODE_RAY_HPP

#include "../MyMath/vec3.hpp"

struct ray {
    vec3 orig;
    vec3 dir;
    double tm = 0;	//time the ray exists at

    ray() = default;
    ray(const vec3& origin, const vec3& direction, const double time = 0.0) : orig(origin), dir(direction), tm(time) {}

    [[nodiscard]] inline vec3 origin() const {return orig;}
    [[nodiscard]] inline vec3 direction() const {return dir;}
    [[nodiscard]] inline double time() const {return tm;}

    [[nodiscard]] vec3 at(const double t) const {
        return orig + t*dir;
    }
};


#endif //CODE_RAY_HPP
