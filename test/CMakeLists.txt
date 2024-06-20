# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2021 Vaibhav Thakkar

# Contributed and/or modified by Vaibhav Thakkar
# Contributed and/or modified by Ioannis Iakovidis
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

enable_testing()

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# Locate Intel MKL root (in case it is enabled)

if (APPLE)
	set(MKLROOT /opt/intel/oneapi/mkl/latest)
elseif(UNIX)
	set(MKLROOT /opt/intel/oneapi/mkl/latest)
endif()


option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../external" OFF)
option(BUILTIN_AUTODIFF "Use autodiff from ../external" ON)
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

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)

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

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

include("../external/cmake-files/Autodiff.cmake")
GetAutodiff()

include("../external/cmake-files/Eigen.cmake")
GetEigen()

include("../external/cmake-files/Boost.cmake")
GetBoost()

include("../external/cmake-files/LPSolve.cmake")
GetLPSolve()

include("../external/cmake-files/QD.cmake")
GetQD()

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O1        # O0 (or no) optimization takes too much time and causes CircleCI test failure.
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

if (BUILTIN_AUTODIFF)
    include_directories (BEFORE ../../external/_deps/Autodiff)
else ()
    include_directories(BEFORE /usr/local/include)
endif(BUILTIN_AUTODIFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

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

include_directories (BEFORE ../external)
include_directories (BEFORE ../include)

#for Eigen
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  add_compile_options(-D "EIGEN_NO_DEBUG")
else ()
  add_compile_definitions("EIGEN_NO_DEBUG")
endif ()

add_definitions(${CMAKE_CXX_FLAGS} "-g")  # enable debuger
#add_definitions(${CMAKE_CXX_FLAGS} "-Wint-in-bool-context")
#add_definitions(${CMAKE_CXX_FLAGS} "-Wall")

add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17") #enable the c++17 support needed by autodiff
#add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
add_definitions(${CMAKE_CXX_FLAGS} "-DMKL_ILP64")
#add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
#add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )



add_library(test_main OBJECT test_main.cpp)





#add_executable (benchmarks_crhmc benchmarks_crhmc.cpp )
#add_executable (benchmarks_crhmc_sampling benchmarks_crhmc_sampling.cpp )



add_executable (test_internal_points test_internal_points.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME test_max_ball
          COMMAND test_internal_points -tc=test_max_ball)
add_test(NAME test_feasibility_point
          COMMAND test_internal_points -tc=test_feasibility_point)
add_test(NAME test_analytic_center
          COMMAND test_internal_points -tc=test_analytic_center)


set(ADDITIONAL_FLAGS "-march=native -DSIMD_LEN=0 -DTIME_KEEPING")

#set_target_properties(benchmarks_crhmc
#                              PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS})
#set_target_properties(benchmarks_crhmc_sampling
#									
TARGET_LINK_LIBRARIES(test_internal_points lp_solve ${MKL_LINK} coverage_config)

#TARGET_LINK_LIBRARIES(benchmarks_crhmc_sampling lp_solve ${MKL_LINK} QD_LIB coverage_config)
#TARGET_LINK_LIBRARIES(benchmarks_crhmc lp_solve ${MKL_LINK} QD_LIB  coverage_config)

