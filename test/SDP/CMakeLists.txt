# VolEsti (volume computation and sampling library)
# Copyright (c) 20012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

enable_testing()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#if(NOT DEFINED BOOST)
#  message(FATAL_ERROR "This program requires the boost library, and will not be compiled.")
#else()
#  message(STATUS "BOOST is DEFINED")
#endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

include_directories (BEFORE ../../external/Eigen)
include_directories (BEFORE ../../external)
include_directories (BEFORE ../../external/minimum_ellipsoid)
#include_directories (BEFORE ../include/cartesian_geom)
#include_directories (BEFORE ../include/convex_bodies)
include_directories (BEFORE ../../external/LPsolve_src/run_headers)
include_directories (BEFORE ../../external/boost)
#include_directories (BEFORE BOOST)
include_directories (BEFORE ../../include/generators)
include_directories (BEFORE ../../include/volume)
include_directories (BEFORE ../../include)
include_directories (BEFORE ../../include/convex_bodies)
include_directories (BEFORE ../../include/convex_bodies/spectrahedra)
include_directories (BEFORE ../../include/convex_bodies/optimization)
include_directories (BEFORE ../../include/convex_bodies/matrices)
include_directories (BEFORE ../../include/annealing)
include_directories (BEFORE ../../include/samplers)
include_directories (BEFORE ../../include/lp_oracles)
include_directories (BEFORE ../../include/misc)

#for Eigen
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  add_compile_options(-D "EIGEN_NO_DEBUG")
else ()
  add_compile_definitions("EIGEN_NO_DEBUG")
endif ()

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
add_definitions(${CMAKE_CXX_FLAGS} "-g")  # enable debuger
#add_definitions(${CMAKE_CXX_FLAGS} "-Wint-in-bool-context")
#add_definitions(${CMAKE_CXX_FLAGS} "-Wall")

add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
#add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
#add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
#add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )


# Find LAPACK and BLAS
# OPENBLAS or ( ( SystemOpenblas or BLAS) and LAPACK)
## prefer local openblas
find_library(OPENBLAS_LIB openblas PATHS external NO_DEFAULT_PATH)
IF (OPENBLAS_LIB)
    set(LAPACK_LIBRARIES ${OPENBLAS_LIB}) #local openblas has lapack build in
    message( STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}" )

    # ARPACK
    find_library(ARPACK_LIB arpack PATHS external NO_DEFAULT_PATH)
    IF (ARPACK_LIB)
        message( STATUS "ARPACK_LIB found: ${ARPACK_LIB}" )

        add_executable (HMC HMC.cpp)


        TARGET_LINK_LIBRARIES(HMC ${LAPACK_LIBRARIES})
        TARGET_LINK_LIBRARIES(HMC ${ARPACK_LIB})

    ELSE()
        message(FATAL_ERROR "This program requires the arpack library, and will not be compiled.")
    ENDIF()

ELSE()
    message(FATAL_ERROR "This program requires the openblas library, and will not be compiled.")
ENDIF()
