#define MOVING_WALL

#include "Fluid_flow/create_flow.hpp"

#include "Examples/calc_derivs.hpp"
#include "Examples/big_vec_derivs.hpp"

const double wx = 3;
const double wy = 4;
const double wz = 5;
const double max_t = 1;
const double Re = 150;

#define NEW_FLOW

int main() {

    calc_derivs_ex();
    big_veg_derivs_ex();



    constexpr output_settings o{};
    std::vector<vec3> pos;
    std::vector<unsigned> inds;
    std::vector<double> mass;
    std::vector<vec3> vels;
    std::vector<vec3> norms;

    constexpr double mass1 = 0.01;

    const vec3 vel_cm = vec3(0.75, 0, 0);
    const vec3 w{};

    const double z_mid = wz/2;
    const double y_mid = wy/2;
    const double x_mid = wx/2;

    //left/right rectangles
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                pos.emplace_back( (2*i-1)*wx/4 + x_mid, (2*j-1)*wy/4 + y_mid, (2*k-1)*wz/4 + z_mid );
                mass.push_back(mass1);
                norms.emplace_back((2*i-1), 0, 0);
            }
        }
    }

    //top/bottom rectangles
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                pos.emplace_back( (2*i-1)*wx/4 + x_mid, (2*j-1)*wy/4 + y_mid, (2*k-1)*wz/4 + z_mid );
                mass.push_back(mass1);
                norms.emplace_back(0, (2*j-1), 0);
            }
        }
    }

    //front/back rectangles
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
                pos.emplace_back( (2*i-1)*wx/4 + x_mid, (2*j-1)*wy/4 + y_mid, (2*k-1)*wz/4 + z_mid );
                mass.push_back(mass1);
                norms.emplace_back(0, 0, (2*k-1));
            }
        }
    }

    for (int i = 0; i < 6; i++) {
        inds.push_back(4*i);    //lower left
        inds.push_back(4*i+1);  //lower right
        inds.push_back(4*i+2);  //upper left

        inds.push_back(4*i+1);  //lower right
        inds.push_back(4*i+2);  //upper left
        inds.push_back(4*i+3);  //upper right
    }

    //mesh m(pos, inds, mass, norms, vec3(0.75,0,0), vec3(0,0,0));
    mesh m(pos, inds, mass, norms, vec3(0.0,0,0), vec3(0,0,0));
    body b(m);
    solve_flow<1000, 128, 128, 128>(&b, o, max_t, Re, wx, wy, wz);


    return 0;
}
