# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2026 Vissarion Fisikopoulos
# Copyright (c) 2018-2026 Apostolos Chalkis
# Copyright (c) 2026      Dimitrios Pavlou

# Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

# Licensed under GNU LGPL.3, see LICENCE file

project (VolEsti)

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

add_definitions(-DDISABLE_NLP_ORACLES)

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/Spectra.cmake")
GetSpectra()

include("../../external/cmake-files/HiGHS.cmake")
GetHighs()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

include("../../external/cmake-files/SuiteSparse.cmake")
GetSuiteSparse()

include("../../external/cmake-files/JSON.cmake")
GetJSON()

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

include_directories (BEFORE ../../external)
include_directories (BEFORE ../../include)
include_directories (BEFORE ../../include/convex_bodies)
include_directories (BEFORE ../../include/misc)
include_directories (BEFORE ../../include/preprocess)

# for Eigen
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
endif ()


add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17")  # enable C++17 standard
add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
add_definitions(${CMAKE_CXX_FLAGS} "-g")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
add_executable(metabolic_simplification metabolic_simplification.cpp)
TARGET_LINK_LIBRARIES(metabolic_simplification highs::highs lp_solve nlohmann_json::nlohmann_json ${SUITESPARSE_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MKL_LINK})