# Licensed under GNU LGPL.3, see LICENCE file
# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2021 Vissarion Fisikopoulos
# Copyright (c) 2018-2021 Apostolos Chalkis

# Contributed and/or modified by Huu Phuoc Le

# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )


CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(MKLROOT /opt/intel/oneapi/mkl/latest)
message("MKLROOT = " ${MKLROOT})
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)
option(USE_MKL "Use MKL library to build eigen" OFF)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
else()
  find_library(IFOPT NAMES libifopt_core.so PATHS /usr/local/lib)
  find_library(IFOPT_IPOPT NAMES libifopt_ipopt.so PATHS /usr/local/lib)
  find_library(GMP NAMES libgmp.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(MPSOLVE NAMES libmps.so PATHS /usr/local/lib)
  find_library(PTHREAD NAMES libpthread.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(FFTW3 NAMES libfftw3.so.3 PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)

  if (NOT IFOPT)

    message(FATAL_ERROR "This program requires the ifopt library, and will not be compiled.")

  elseif (NOT GMP)

    message(FATAL_ERROR "This program requires the gmp library, and will not be compiled.")

  elseif (NOT MPSOLVE)

    message(FATAL_ERROR "This program requires the mpsolve library, and will not be compiled.")

  elseif (NOT FFTW3)

    message(FATAL_ERROR "This program requires the fftw3 library, and will not be compiled.")

  else()
    message(STATUS "Library ifopt found: ${IFOPT}")
    message(STATUS "Library gmp found: ${GMP}")
    message(STATUS "Library mpsolve found: ${MPSOLVE}")
    message(STATUS "Library fftw3 found:" ${FFTW3})

  endif(NOT IFOPT)

endif(DISABLE_NLP_ORACLES)

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

if (BUILTIN_EIGEN)
    include_directories (BEFORE ../../external/_deps/Eigen)
else ()
    include_directories(BEFORE /usr/include/eigen3)
endif(BUILTIN_EIGEN)

if (USE_MKL)
  find_library(BLAS NAME libblas.so
  PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu /usr/lib64/atlas /usr/lib/atlas /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib
  )
  message(STATUS "BLAS library found: ${BLAS}")

  message(STATUS "MKL library found: ${MKLROOT}")
  include_directories (BEFORE ${MKLROOT}/include)
  set(PROJECT_LIBS ${BLAS_LIBRARIES})
  add_definitions(-DEIGEN_USE_MKL_ALL)

  # Set linking for MKL : may affect performance
  set(MKL_LINK "-L${MKLROOT}/lib/intel64 -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lgomp  -lpthread -llapack -lblas")
endif(USE_MKL)

# Find lpsolve library
find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")
else ()
  message(STATUS "lp_solve library found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../../external)
  include_directories (BEFORE ../../external/minimum_ellipsoid)
  include_directories (BEFORE ../../include/generators)
  include_directories (BEFORE ../../include/volume)
  include_directories (BEFORE ../../include)
  include_directories (BEFORE ../../include/lp_oracles)
  include_directories (BEFORE ../../include/nlp_oracles)
  include_directories (BEFORE ../../include/convex_bodies)
  include_directories (BEFORE ../../include/random_walks)
  include_directories (BEFORE ../../include/annealing)
  include_directories (BEFORE ../../include/ode_solvers)
  include_directories (BEFORE ../../include/root_finders)
  include_directories (BEFORE ../../include/samplers)
  include_directories (BEFORE ../../include/misc)
  include_directories (BEFORE ../../include/optimization)

  # for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()

  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CMAKE_CXX_FLAGS} "-DMKL_ILP64")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")

  add_executable (sampler sampler.cpp)
  TARGET_LINK_LIBRARIES(sampler ${LP_SOLVE} ${BLAS} ${MKL_LINK})

endif()
