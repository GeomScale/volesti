# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
# Contributed and/or modified by Angelos Korakitis, as part of Google Summer of 2025 program.
# Licensed under GNU LGPL.3, see LICENCE file

# project( VolEsti )

# enable_testing()

# cmake_minimum_required(VERSION 4.0)

# set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# if(COMMAND cmake_policy)
#        cmake_policy(SET CMP0003 NEW)
# endif(COMMAND cmake_policy)

# #if(NOT DEFINED BOOST)
# #  message(FATAL_ERROR "This program requires the boost library, and will not be compiled.")
# #else()
# #  message(STATUS "BOOST is DEFINED")
# #endif()


# set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

# include_directories (BEFORE ../../external/Eigen)
# include_directories (BEFORE ../../external)
# include_directories (BEFORE ../../external/minimum_ellipsoid)
# #include_directories (BEFORE ../include/cartesian_geom)
# #include_directories (BEFORE ../include/convex_bodies)
# include_directories (BEFORE ../../external/LPsolve_src/run_headers)
# include_directories (BEFORE ../../external/boost)
# #include_directories (BEFORE BOOST)
# include_directories (BEFORE ../../include/generators)
# include_directories (BEFORE ../../include/volume)
# include_directories (BEFORE ../../include)
# include_directories (BEFORE ../../include/convex_bodies)
# include_directories (BEFORE ../../include/convex_bodies/spectrahedra)
# include_directories (BEFORE ../../include//optimization)
# include_directories (BEFORE ../../include/matrix_operations)
# include_directories (BEFORE ../../include/annealing)
# include_directories (BEFORE ../../include/samplers)
# include_directories (BEFORE ../../include/lp_oracles)
# include_directories (BEFORE ../../include/misc)

# #for Eigen
# if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
#   add_compile_options(-D "EIGEN_NO_DEBUG")
# else ()
#   add_compile_definitions("EIGEN_NO_DEBUG")
# endif ()

# add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
# add_definitions(${CMAKE_CXX_FLAGS} "-g")  # enable debuger
# #add_definitions(${CMAKE_CXX_FLAGS} "-Wint-in-bool-context")
# #add_definitions(${CMAKE_CXX_FLAGS} "-Wall")

# add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
# #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
# add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
# add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
# add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
# #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
# #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )


# # Find LAPACK and BLAS
# # OPENBLAS or ( ( SystemOpenblas or BLAS) and LAPACK)
# ## prefer local openblas
# find_library(OPENBLAS_LIB openblas PATHS external NO_DEFAULT_PATH)
# IF (OPENBLAS_LIB)
#     set(LAPACK_LIBRARIES ${OPENBLAS_LIB}) #local openblas has lapack build in

#     find_package( Threads REQUIRED )
#     set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#     message( STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}" )


#     # ARPACK
#     find_library(ARPACK_LIB arpack PATHS external NO_DEFAULT_PATH)
#     IF (ARPACK_LIB)
#         message( STATUS "ARPACK_LIB found: ${ARPACK_LIB}" )

#         # Query gfortran to get the libgfortran path
#        FIND_LIBRARY(GFORTRAN_LIB NAMES libgfortran.so PATHS /usr/lib/gcc/x86_64-linux-gnu/8/)

#         IF (CMAKE_Fortran_COMPILER)
#             EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.a
#                     OUTPUT_VARIABLE _libgfortran_path
#                     OUTPUT_STRIP_TRAILING_WHITESPACE
#                     )
#             IF (NOT EXISTS ${_libgfortran_path})
#                 EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.so
#                         OUTPUT_VARIABLE _libgfortran_path
#                         OUTPUT_STRIP_TRAILING_WHITESPACE
#                         )
#             ENDIF ()
#         ENDIF()

#         IF(EXISTS ${_libgfortran_path})
#             get_filename_component(GFORTRAN_PATH ${_libgfortran_path} PATH)
#             find_library(GFORTRAN_LIB gfortran PATHS ${GFORTRAN_PATH})
#         ELSE()
#             # if libgfortran wasn't found at this point, the installation is probably broken
#             # Let's try to find the library nonetheless.
#             FIND_LIBRARY(GFORTRAN_LIB gfortran)
#         ENDIF()

#         IF (GFORTRAN_LIB)
#             message( STATUS "GFORTRAN_LIB found: ${GFORTRAN_LIB}" )

#             add_library(test_main OBJECT test_main.cpp)

#             add_executable(test_solve_sdp test_solve_sdp.cpp $<TARGET_OBJECTS:test_main>)
#             add_test(NAME test_solve_sdp COMMAND test_solve_sdp -tc=test_solve_sdp)
#             TARGET_LINK_LIBRARIES(test_solve_sdp ${ARPACK_LIB} ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})

#             add_executable(test_spec_oracles test_spec_oracles.cpp $<TARGET_OBJECTS:test_main>)
#             add_test(NAME test_spec_oracles COMMAND test_spec_oracles -tc=test_spec_oracles)
#             TARGET_LINK_LIBRARIES(test_spec_oracles ${ARPACK_LIB} ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})

#         ELSE()
#             MESSAGE(STATUS "gfortran is required but it could not be found")
#         ENDIF ()


#     ELSE()
#         message(FATAL_ERROR "This program requires the arpack library, and will not be compiled.")
#     ENDIF()

# ELSE()
#     message(FATAL_ERROR "This program requires the openblas library, and will not be compiled.")
# ENDIF()


cmake_minimum_required(VERSION 3.15)

project(VolEsti LANGUAGES CXX Fortran)

# --------------------------------------------------
#  Options users can override from the command line
# --------------------------------------------------
option(USE_SYSTEM_EIGEN "Prefer system‑installed Eigen3 if available" ON)
set(EIGEN3_ROOT "" CACHE PATH "Path to Eigen3 *parent of Eigen/* (e.g. /opt/homebrew/opt/eigen/include/eigen3)")
set(OPENBLAS_ROOT "" CACHE PATH "Preferred OpenBLAS prefix (parent of lib/libopenblas.*)")
set(ARPACK_ROOT "" CACHE PATH "Preferred ARPACK prefix (parent of lib/libarpack.*)")

# --------------------------------------------------
#  Build‑type & compile flags
# --------------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_definitions(EIGEN_NO_DEBUG)
add_compile_options($<$<CONFIG:Release>:-O3>)
add_compile_options(-DBOOST_NO_AUTO_PTR)

# --------------------------------------------------
#  Dependencies
# --------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH ${OPENBLAS_ROOT} ${ARPACK_ROOT} ${EIGEN3_ROOT})

# --- Eigen3 (header‑only) ----------------------------------------------------
if (USE_SYSTEM_EIGEN)
  find_package(Eigen3 QUIET)
endif()

if (NOT Eigen3_FOUND)
  if (EXISTS ${CMAKE_SOURCE_DIR}/external/eigen3/Eigen)
    set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/eigen3)
  elseif(EIGEN3_ROOT AND EXISTS ${EIGEN3_ROOT}/Eigen)
    set(EIGEN3_INCLUDE_DIR ${EIGEN3_ROOT})
  else()
    message(FATAL_ERROR "Eigen3 headers not found. Either install Eigen3 and set USE_SYSTEM_EIGEN=ON or place a copy under external/eigen3/.")
  endif()
else()
  set(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIRS})
endif()
include_directories(BEFORE ${EIGEN3_INCLUDE_DIR})

# --- OpenBLAS  --------------------------------------------------------------
find_library(OPENBLAS_LIB openblas
             PATHS ${OPENBLAS_ROOT} ${CMAKE_SOURCE_DIR}/external
             PATH_SUFFIXES lib
             NO_DEFAULT_PATH)
if (NOT OPENBLAS_LIB)
  find_library(OPENBLAS_LIB openblas) # system paths
endif()
if (NOT OPENBLAS_LIB)
  message(FATAL_ERROR "OpenBLAS library not found. Set OPENBLAS_ROOT or install OpenBLAS.")
endif()
set(LAPACK_LIBRARIES ${OPENBLAS_LIB})
find_package(Threads REQUIRED)
list(APPEND LAPACK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# --- ARPACK -----------------------------------------------------------------
find_library(ARPACK_LIB arpack
             PATHS ${ARPACK_ROOT} ${CMAKE_SOURCE_DIR}/external
             PATH_SUFFIXES lib
             NO_DEFAULT_PATH)
if (NOT ARPACK_LIB)
  find_library(ARPACK_LIB arpack) # system paths
endif()
if (NOT ARPACK_LIB)
  message(FATAL_ERROR "ARPACK library not found. Set ARPACK_ROOT or install ARPACK.")
endif()

# --- libgfortran ------------------------------------------------------------
if (CMAKE_Fortran_COMPILER)
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.dylib
                  OUTPUT_VARIABLE _libgfortran_path
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (_libgfortran_path AND EXISTS ${_libgfortran_path})
    get_filename_component(GFORTRAN_LIB ${_libgfortran_path} ABSOLUTE)
  else()
    find_library(GFORTRAN_LIB gfortran)
  endif()
else()
  find_library(GFORTRAN_LIB gfortran)
endif()
if (NOT GFORTRAN_LIB)
  message(FATAL_ERROR "libgfortran not found – install GCC with gfortran.")
endif()

# --- Boost (headers only for now) ------------------------------------------
# Users may point to a custom Boost via BOOST_ROOT or standard CMAKE_PREFIX_PATH.
set(BOOST_ROOT "" CACHE PATH "Prefix where Boost is installed (parent of include/boost)")
list(APPEND CMAKE_PREFIX_PATH ${BOOST_ROOT})
find_package(Boost 1.65 REQUIRED)
include_directories(BEFORE ${Boost_INCLUDE_DIRS})

# --------------------------------------------------
#  Include directories for project headers
# --------------------------------------------------
# If this CMakeLists.txt is run from test/SDP/ the real project root is two
# levels up (../..). Compute that so header paths work no matter where we
# configure from.
get_filename_component(PROJECT_TOP ${CMAKE_SOURCE_DIR}/../.. ABSOLUTE)

include_directories(BEFORE
    ${PROJECT_TOP}/external
    ${PROJECT_TOP}/external/minimum_ellipsoid
    ${PROJECT_TOP}/external/LPsolve_src/run_headers
    ${PROJECT_TOP}/external/boost
    ${PROJECT_TOP}/include
    ${PROJECT_TOP}/include/cartesian_geom
    ${PROJECT_TOP}/include/generators
    ${PROJECT_TOP}/include/volume
    ${PROJECT_TOP}/include/convex_bodies
    ${PROJECT_TOP}/include/convex_bodies/spectrahedra
    ${PROJECT_TOP}/include/optimization
    ${PROJECT_TOP}/include/matrix_operations
    ${PROJECT_TOP}/include/annealing
    ${PROJECT_TOP}/include/samplers
    ${PROJECT_TOP}/include/lp_oracles
    ${PROJECT_TOP}/include/misc)

# --------------------------------------------------
#  Unit tests in test/SDP
# --------------------------------------------------
add_library(test_main OBJECT test_main.cpp)

function(volesti_add_test name)
  add_executable(${name} ${ARGV1} $<TARGET_OBJECTS:test_main>)
  add_test(NAME ${name} COMMAND ${name} -tc=${name})
  target_link_libraries(${name} PRIVATE ${ARPACK_LIB} ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})
endfunction()

volesti_add_test(test_solve_sdp test_solve_sdp.cpp)
volesti_add_test(test_solve_sparse_sdp test_solve_sparse_sdp.cpp)
volesti_add_test(test_spec_oracles test_spec_oracles.cpp)
