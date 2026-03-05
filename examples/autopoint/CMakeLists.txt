# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis

project( VolEsti )

CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../external" OFF)
option(USE_MKL "Use MKL library to build eigen" OFF)

# NEW: Optional HiGHS solver
option(USE_HIGHS "Enable HiGHS LP solver" ON)

set(MKLROOT /opt/intel/oneapi/mkl/2023.0.0)
message("MKLROOT ${MKLROOT}")

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif()

############################################
# NLP ORACLES
############################################

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
else()
  find_library(IFOPT NAMES libifopt_core.so PATHS /usr/local/lib)
  find_library(IFOPT_IPOPT NAMES libifopt_ipopt.so PATHS /usr/local/lib)
  find_library(GMP NAMES libgmp.so PATHS /usr/lib/x86_64-linux-gnu)
  find_library(MPSOLVE NAMES libmps.so PATHS /usr/local/lib)
  find_library(PTHREAD NAMES libpthread.so PATHS /usr/lib/x86_64-linux-gnu)
  find_library(FFTW3 NAMES libfftw3.so.3 PATHS /usr/lib/x86_64-linux-gnu)

  if (NOT IFOPT)
      message(FATAL_ERROR "This program requires the ifopt library.")
  endif()
endif()

############################################
# External dependencies
############################################

option(BUILTIN_AUTODIFF "Use autodiff from ../external" OFF)

include("../../external/cmake-files/Autodiff.cmake")
GetAutodiff()

if (BUILTIN_AUTODIFF)
    include_directories(BEFORE ../../external/_deps/Autodiff)
else ()
    include_directories(BEFORE /usr/local/include)
endif()

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

############################################
# REMOVE HARD LPSOLVE DEPENDENCY
############################################

# OLD (removed)
# include("../../external/cmake-files/LPSolve.cmake")
# GetLPSolve()

# OLD FAIL BLOCK REMOVED
# find_library(LP_SOLVE ...)
# if (NOT LP_SOLVE) message(FATAL_ERROR)

############################################
# OPTIONAL HIGHS SOLVER
############################################

if(USE_HIGHS)
    find_package(highs QUIET)

    if(highs_FOUND)
        message(STATUS "HiGHS solver found")
        set(LP_SOLVER_LIB highs::highs)
        add_definitions(-DVOLESTI_USE_HIGHS)
    else()
        message(WARNING "HiGHS not found. LP functionality may be limited.")
    endif()
endif()

############################################
# INCLUDE DIRECTORIES
############################################

include_directories(BEFORE ../../external)
include_directories(BEFORE ../../external/minimum_ellipsoid)
include_directories(BEFORE ../../include)
include_directories(BEFORE ../../include/generators)
include_directories(BEFORE ../../include/volume)
include_directories(BEFORE ../../include/lp_oracles)
include_directories(BEFORE ../../include/nlp_oracles)
include_directories(BEFORE ../../include/convex_bodies)
include_directories(BEFORE ../../include/random_walks)
include_directories(BEFORE ../../include/annealing)
include_directories(BEFORE ../../include/ode_solvers)
include_directories(BEFORE ../../include/root_finders)
include_directories(BEFORE ../../include/samplers)
include_directories(BEFORE ../../include/misc)
include_directories(BEFORE ../../include/optimization)

############################################
# COMPILER FLAGS
############################################

add_definitions("-std=c++17")
add_definitions("-O3")
add_definitions("-pthread")
add_definitions("-lm")
add_definitions("-ldl")
add_definitions("-DBOOST_NO_AUTO_PTR")
add_definitions("-D_REENTRANT")

############################################
# EXECUTABLES
############################################

add_executable(userDefinedFunction_manual_gradient userDefinedFunction_manual_gradient.cpp)
target_link_libraries(userDefinedFunction_manual_gradient ${LP_SOLVER_LIB})

add_executable(userDefinedFunction_autopoint userDefinedFunction_autopoint.cpp)
target_link_libraries(userDefinedFunction_autopoint ${LP_SOLVER_LIB})

add_executable(Gaussian_mixture_autopoint Gaussian_mixture_autopoint.cpp readData.cpp)
target_link_libraries(Gaussian_mixture_autopoint ${LP_SOLVER_LIB})

add_executable(MultidimensionalGaussian_mixture_autopoint MultidimensionalGaussian_mixture_autopoint.cpp readData.cpp)
target_link_libraries(MultidimensionalGaussian_mixture_autopoint ${LP_SOLVER_LIB})