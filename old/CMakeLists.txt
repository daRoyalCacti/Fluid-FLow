cmake_minimum_required(VERSION 3.19)
project(Code)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

add_library(MyMath MyMath/vec3.cpp MyMath/vec3.hpp MyMath/ops.cpp MyMath/ops.hpp)
add_library(Rigid_body Rigid_body/body.cpp Rigid_body/body.hpp)
add_executable(Code main.cpp)


find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Code PUBLIC OpenMP::OpenMP_CXX)
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Code PUBLIC Eigen3::Eigen)
target_link_libraries(MyMath PUBLIC Eigen3::Eigen)
target_link_libraries(Rigid_body PUBLIC Eigen3::Eigen)

target_link_libraries(Code PUBLIC MyMath)
target_link_libraries(Code PUBLIC Rigid_body)

target_link_libraries(Rigid_body PUBLIC MyMath)


set(warning_flags -Werror -Wpedantic -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -Wmissing-include-dirs -Wunknown-pragmas -Wduplicated-cond -Wfloat-equal -Wshadow -Wunsafe-loop-optimizations -Wpacked -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wmissing-noreturn -Wsuggest-attribute=malloc -Wsuggest-attribute=format -Wmissing-format-attribute -Wsuggest-attribute=cold>)
add_compile_options("$<$<CONFIG:Debug>: ${warning_flags} ")
add_compile_options("$<$<CONFIG:Release>: ${warning_flags} -Ofast -fno-math-errno -funsafe-math-optimizations -ffinite-math-only")

