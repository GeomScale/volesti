# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2025 Vissarion Fisikopoulos
# Copyright (c) 2018-2025 Apostolos Chalkis
# Copyright (c) 2020-2025 Elias Tsigaridas

# Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
# Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

# Licensed under GNU LGPL.3, see LICENCE file

cmake_minimum_required(VERSION 3.15)

project(VolEsti LANGUAGES CXX Fortran)

include(FetchContent)

# --------------------------------------------------
#  Options users can override from the command line
# --------------------------------------------------
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
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_definitions(EIGEN_NO_DEBUG)
add_compile_options($<$<CONFIG:Release>:-O3>)
add_compile_options(-DBOOST_NO_AUTO_PTR)

# --------------------------------------------------
#  Compute project root directory first
# --------------------------------------------------
get_filename_component(PROJECT_TOP ${CMAKE_SOURCE_DIR}/../.. ABSOLUTE)

# --------------------------------------------------
#  Dependencies
# --------------------------------------------------
if (APPLE)
  list(APPEND CMAKE_PREFIX_PATH 
    /opt/homebrew/opt/openblas 
    /opt/homebrew/opt/arpack-ng
    /opt/homebrew
    /usr/local/opt/openblas
    /usr/local/opt/arpack-ng
    /usr/local
  )
endif()

list(APPEND CMAKE_PREFIX_PATH ${OPENBLAS_ROOT} ${ARPACK_ROOT})

# --- Eigen3 (header‑only) ----------------------------------------------------
include(${PROJECT_TOP}/external/cmake-files/Eigen.cmake)
GetEigen()

# --- Spectra (header‑only) ---------------------------------------------------
include(${PROJECT_TOP}/external/cmake-files/Spectra.cmake)
GetSpectra()

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
set(BOOST_ROOT "" CACHE PATH "Prefix where Boost is installed (parent of include/boost)")
list(APPEND CMAKE_PREFIX_PATH ${BOOST_ROOT})
find_package(Boost 1.65 REQUIRED)
include_directories(BEFORE ${Boost_INCLUDE_DIRS})

# --------------------------------------------------
#  Include directories for project headers
# --------------------------------------------------
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
    ${PROJECT_TOP}/include/misc
    ${PROJECT_TOP}/include/preprocess/spectrahedron)

# --------------------------------------------------
#  Unit tests
# --------------------------------------------------
add_library(test_main OBJECT test_main.cpp)

function(volesti_add_test name)
  add_executable(${name} ${ARGV1} $<TARGET_OBJECTS:test_main>)
  add_test(NAME ${name} COMMAND ${name} -tc=${name})
  target_link_libraries(${name} PRIVATE ${ARPACK_LIB} ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})
endfunction()

volesti_add_test(test_solve_sdp test_solve_sdp.cpp)
volesti_add_test(test_spec_oracles test_spec_oracles.cpp)
volesti_add_test(test_solve_sparse_sdp test_solve_sparse_sdp.cpp)