# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2024 Vissarion Fisikopoulos
# Copyright (c) 2018-2024 Apostolos Chalkis
#
# CMakeLists.txt for building Octave MEX interfaces
#
# Licensed under GNU LGPL.3, see LICENSE file

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

# Find Octave
find_program(OCTAVE_CONFIG octave-config PATHS
  /usr/bin
  /usr/local/bin
  /opt/local/bin
  ENV PATH
)

if(NOT OCTAVE_CONFIG)
  message(FATAL_ERROR "octave-config not found. Please install Octave and ensure it's in your PATH.")
endif()

# Get Octave paths
execute_process(
  COMMAND ${OCTAVE_CONFIG} -p OCTINCLUDEDIR
  OUTPUT_VARIABLE OCTAVE_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${OCTAVE_CONFIG} -p OCTLIBDIR
  OUTPUT_VARIABLE OCTAVE_LIB_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${OCTAVE_CONFIG} -p OCTFILESDIR
  OUTPUT_VARIABLE OCTAVE_OCTFILES_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Octave include directory: ${OCTAVE_INCLUDE_DIR}")
message(STATUS "Octave lib directory: ${OCTAVE_LIB_DIR}")
message(STATUS "Octave octfiles directory: ${OCTAVE_OCTFILES_DIR}")

# Include directories
include_directories(${OCTAVE_INCLUDE_DIR})
include_directories(../../external)
include_directories(../../external/minimum_ellipsoid)
include_directories(../../include)
include_directories(../../include/cartesian_geom)
include_directories(../../include/convex_bodies)
include_directories(../../include/generators)
include_directories(../../include/volume)
include_directories(../../include/sampling)
include_directories(../../include/random_walks)
include_directories(../../include/lp_oracles)
include_directories(../../include/preprocess)
include_directories(../../include/misc)

# Eigen
option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)
include("../../external/cmake-files/Eigen.cmake")
GetEigen()

# Boost
include("../../external/cmake-files/Boost.cmake")
GetBoost()

# lp_solve
include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

# Compiler flags
add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17")
add_definitions(${CMAKE_CXX_FLAGS} "-O3")
add_definitions(${CMAKE_CXX_FLAGS} "-DBOOST_NO_AUTO_PTR")
add_definitions(${CMAKE_CXX_FLAGS} "-DEIGEN_NO_DEBUG")

if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  add_compile_options(-D "EIGEN_NO_DEBUG")
else ()
  add_compile_definitions("EIGEN_NO_DEBUG")
endif ()

# Source files
set(MEX_SOURCES
  volesti_volume_mex.cpp
  volesti_sample_mex.cpp
)

# Build MEX files
foreach(MEX_FILE ${MEX_SOURCES})
  get_filename_component(MEX_NAME ${MEX_FILE} NAME_WE)
  
  add_library(${MEX_NAME} SHARED ${MEX_FILE})
  target_link_libraries(${MEX_NAME} PUBLIC lp_solve)
  
  # Set output directory
  set_target_properties(${MEX_NAME} PROPERTIES
    PREFIX ""
    SUFFIX ".oct"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../inst"
  )
endforeach()

