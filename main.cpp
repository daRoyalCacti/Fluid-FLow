#define UPDATE_VECS_CHECK_RESULTS_LOG

#include "Fluid_flow/create_flow.hpp"

//#include "Examples/calc_derivs.hpp"
//#include "Examples/big_vec_derivs.hpp"
//#include "Examples/interp.hpp"

#define HORIZONTAL
#define ROTATE

#ifdef HORIZONTAL
constexpr double wx = 7;
constexpr double wy = 6;
constexpr double wz = 5;
#else
constexpr double wx = 6;
constexpr double wy = 7;
constexpr double wz = 5;
#endif
constexpr double max_t = 1;
constexpr double Re = 7069; //150;  //https://www.grc.nasa.gov/WWW/k-12/airplane/reynolds.html




int main() {
    //big_veg_derivs_ex();
    //calc_derivs_ex();
    //interp_ex();

//*
    constexpr output_settings o{};
    std::vector<vec3> pos;
    std::vector<unsigned> inds;
    std::vector<double> mass;
    std::vector<vec3> vels;
    std::vector<vec3> norms;

    constexpr double mass1 = 0.01;
#ifdef HORIZONTAL
    constexpr vec3 vel_cm = vec3(0.75, 0, 0);
#else
    constexpr vec3 vel_cm = vec3(0, 0.75, 0);//vec3(0.1, 0, 0);//vec3(0, 0, 0);
#endif
#ifdef ROTATE
    constexpr vec3 w = vec3(0,0,1);
#else
    constexpr vec3 w = vec3(0,0,0);
#endif

    const double z_mid = wz/2;
    const double y_mid = wy/2;
    const double x_mid = wx/2;


    constexpr double box_x = 1;
    constexpr double box_y = 1;
    constexpr double box_z = 1;

    //left/right rectangles
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                pos.emplace_back( (2*i-1)*box_x/2 + x_mid, (2*j-1)*box_y/2 + y_mid, (2*k-1)*box_z/2 + z_mid );
                mass.push_back(mass1);
                norms.emplace_back((2*i-1), 0, 0);
            }
        }
    }

    //top/bottom rectangles
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                pos.emplace_back( (2*i-1)*box_x/2 + x_mid, (2*j-1)*box_y/2 + y_mid, (2*k-1)*box_z/2 + z_mid );
                mass.push_back(mass1);
                norms.emplace_back(0, (2*j-1), 0);
            }
        }
    }

    //front/back rectangles
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
                pos.emplace_back( (2*i-1)*box_x/2 + x_mid, (2*j-1)*box_y/2 + y_mid, (2*k-1)*box_z/2 + z_mid );
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
    mesh m(pos, std::move(inds), std::move(mass), std::move(norms), vel_cm, w);
    body b(m);
#ifdef NDEBUG
    solve_flow<1000, 128, 128, 128>(&b, o, max_t, Re, wx, wy, wz);
#else
    solve_flow<1000, 30, 30, 30>(&b, o, max_t, Re, wx, wy, wz);
#endif
//*/

    return 0;
}
