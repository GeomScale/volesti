# volesti benchmarking: Convex Hull Approximation
#
# Demonstrates: inner approximation via convex hull, support function gap
# measurement, convergence rate validation, and body type auto-detection
# from vertex count growth.

project(geometric_inference)
cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-DDISABLE_NLP_ORACLES)

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

include_directories(
    BEFORE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../external
    ${CMAKE_CURRENT_SOURCE_DIR}/../../external/minimum_ellipsoid
)

add_executable(geometric_inference_poc geometric_inference_poc.cpp)
target_link_libraries(geometric_inference_poc PUBLIC lp_solve)
