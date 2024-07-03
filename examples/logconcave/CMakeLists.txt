# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.
# Licensed under GNU LGPL.3, see LICENCE file


project( VolEsti )


CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(MKLROOT /opt/intel/oneapi/mkl/latest)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../external" OFF)
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

find_library(BLAS NAMES libblas.so libblas.dylib PATHS /usr/local/Cellar/lapack/3.9.1_1/lib /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu /usr/local/Cellar/openblas/0.3.15_1/lib /usr/lib)

if (USE_MKL)
  find_library(BLAS NAMES libblas.so libblas.dylib PATHS /usr/local/Cellar/lapack/3.9.1_1/lib /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu /usr/local/Cellar/openblas/0.3.15_1/lib /usr/lib)
  find_library(GFORTRAN NAME libgfortran.dylib PATHS /usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10)
  find_library(LAPACK NAME liblapack.dylib PATHS /usr/lib)
  find_library(OPENMP NAME libiomp5.dylib PATHS /opt/intel/oneapi/compiler/2021.1.1/mac/compiler/lib)

  include_directories (BEFORE ${MKLROOT}/include)
  set(PROJECT_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GFORTRAN_LIBRARIES})
  set(MKL_LINK "-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
  add_definitions(-DEIGEN_USE_MKL_ALL)
else()
  set(MKL_LINK "")
endif(USE_MKL)

# Find lpsolve library
find_library(LP_SOLVE NAMES liblpsolve55.so liblpsolve55.dylib PATHS /usr/lib/lp_solve /usr/local/lib)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")
else ()
  message(STATUS "Library lp_solve found: ${LP_SOLVE} ${BLAS}")

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


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  add_definitions(${CMAKE_CXX_FLAGS} "-DMKL_ILP64")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )

  add_executable (simple_hmc simple_hmc.cpp)
  add_executable (simple_ode simple_ode.cpp)
  add_executable (simple_uld simple_uld.cpp)
  add_executable (simple_sde simple_sde.cpp)
  add_executable (high_dimensional_hmc high_dimensional_hmc.cpp)
  add_executable (high_dimensional_hmc_rand_poly high_dimensional_hmc_rand_poly.cpp)
  add_executable (exponential_exact_hmc exponential_exact_hmc.cpp)
  add_executable (gaussian_exact_hmc gaussian_exact_hmc.cpp)

  #add_executable (adaptive_step_size adaptive_step_size.cpp)

  TARGET_LINK_LIBRARIES(high_dimensional_hmc_rand_poly ${LP_SOLVE} ${BLAS} ${MKL_LINK} )
  TARGET_LINK_LIBRARIES(high_dimensional_hmc ${LP_SOLVE} ${BLAS} ${MKL_LINK}  )
  #TARGET_LINK_LIBRARIES(adaptive_step_size ${LP_SOLVE} ${BLAS}  ${MKL_LINK}  )
  TARGET_LINK_LIBRARIES(simple_ode ${LP_SOLVE} ${BLAS} ${MKL_LINK}  )
  TARGET_LINK_LIBRARIES(simple_sde ${LP_SOLVE} ${BLAS} ${MKL_LINK}  )
  TARGET_LINK_LIBRARIES(simple_uld ${LP_SOLVE} ${BLAS} ${MKL_LINK}   )
  TARGET_LINK_LIBRARIES(simple_hmc ${LP_SOLVE} ${BLAS} ${MKL_LINK} )
  TARGET_LINK_LIBRARIES(exponential_exact_hmc ${LP_SOLVE} ${BLAS} ${MKL_LINK}  )
  TARGET_LINK_LIBRARIES(gaussian_exact_hmc ${LP_SOLVE} ${BLAS} ${MKL_LINK}  )


endif()
