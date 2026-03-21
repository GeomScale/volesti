# VolEsti (volume computation and sampling library)
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
endif(DISABLE_NLP_ORACLES)

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

# Find lpsolve library
find_library(LP_SOLVE NAMES liblpsolve55.so liblpsolve55.a PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")
else ()
  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../../external)
  include_directories (BEFORE ../../external/minimum_ellipsoid)
  include_directories (BEFORE ../../include/generators)
  include_directories (BEFORE ../../include/volume)
  include_directories (BEFORE ../../include)
  include_directories (BEFORE ../../include/convex_bodies)
  include_directories (BEFORE ../../include/random_walks)
  include_directories (BEFORE ../../include/annealing)
  include_directories (BEFORE ../../include/ode_solvers)
  include_directories (BEFORE ../../include/root_finders)
  include_directories (BEFORE ../../include/samplers)
  include_directories (BEFORE ../../include/lp_oracles)
  include_directories (BEFORE ../../include/nlp_oracles)
  include_directories (BEFORE ../../include/misc)
  include_directories (BEFORE ../../include/optimization)
  include_directories (BEFORE ../../include/preprocess)

  # for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()

  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17")  # enable C++17 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")

  add_executable (count_le_orderpolytope count_le_orderpolytope.cpp)
  TARGET_LINK_LIBRARIES(count_le_orderpolytope ${LP_SOLVE})

endif()
