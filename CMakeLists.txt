cmake_minimum_required(VERSION 3.19)
project(Code LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")



add_executable(Code main.cpp MyMath/vec3.hpp MyMath/ops.hpp Rigid_body/body.hpp Examples/Rigid_body_ex.hpp MyMath/finite_difference.hpp Examples/derivatives.hpp MyMath/big_vec.hpp Examples/big_vec_derivs.h MyMath/calc.hpp Examples/calc_derivs.hpp MyMath/big_matrix.hpp Fluid_flow/make_vecs.hpp Fluid_flow/make_mats.hpp Fluid_flow/create_flow.hpp Fluid_flow/solver.hpp MyMath/boundary.hpp Fluid_flow/flow_env.hpp Fluid_flow/timing.hpp Rigid_body/mesh.hpp Rigid_body/triangle.hpp Rigid_body/ray.hpp Fluid_flow/boundary_conditions.hpp Rigid_body/triangle_mesh.hpp Fluid_flow/update_mesh.hpp)


set_target_properties(Code PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Code PUBLIC OpenMP::OpenMP_CXX)
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Code PUBLIC Eigen3::Eigen)

find_package(OpenCL REQUIRED)
target_link_libraries(Code PUBLIC OpenCL::OpenCL)


set(warning_flags -Werror -Wpedantic -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -Wmissing-include-dirs -Wunknown-pragmas -Wduplicated-cond -Wfloat-equal -Wshadow -Wunsafe-loop-optimizations -Wpacked -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wmissing-noreturn -Wsuggest-attribute=malloc -Wsuggest-attribute=format -Wmissing-format-attribute -Wsuggest-attribute=cold>)
add_compile_options("$<$<CONFIG:Debug>: ${warning_flags}")
add_compile_options("$<$<CONFIG:Release>: ${warning_flags} -Ofast -fno-math-errno -funsafe-math-optimizations -ffinite-math-only -march=native -mfma")

add_definitions(-DVIENNACL_WITH_OPENCL -DEIGEN_NO_AUTOMATIC_RESIZING)


