# Licensed under GNU LGPL.3, see LICENCE file
# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2021 Vissarion Fisikopoulos
# Copyright (c) 2018-2021 Apostolos Chalkis
# Contributed and/or modified by Vladimir Necula, as part of Google Summer of Code 2025 program.

# Licensed under GNU LGPL.3, see LICENCE file

cmake_minimum_required(VERSION 3.10)
project(SparseBilliardWalkBenchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
endif()

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

include_directories(${CMAKE_SOURCE_DIR}/../../include)
include_directories(${CMAKE_SOURCE_DIR}/../../external)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(sparse_billiard_benchmark sparse_billiard_walk_benchmark.cpp)

target_link_libraries(sparse_billiard_benchmark 
    lp_solve
    ${Boost_LIBRARIES}
)

target_compile_options(sparse_billiard_benchmark PRIVATE 
    -O3 
    -DEIGEN_NO_DEBUG
    -march=native
    -funroll-loops
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sparse_billiard_benchmark OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(sparse_billiard_benchmark PRIVATE DISABLE_NLP_ORACLES) 