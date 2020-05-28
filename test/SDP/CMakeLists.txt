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
ELSE()
    find_library(OPENBLAS_LIB openblas PATH_SUFFIXES openblas-base)
    IF (OPENBLAS_LIB)
        find_library(LAPACK_LIBRARIES lapack)  #only add liblapack
        IF (NOT LAPACK_LIBRARIES)
            message( STATUS "LAPACK_LIBRARIES not found" )
            SET(ABORT_CONFIG TRUE)
        ENDIF()
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${OPENBLAS_LIB})
    ELSE()
        message( STATUS "OPENBLAS_LIB not found, trying LAPACK/BLAS..." )
        find_package(LAPACK REQUIRED)  # this will also find and add libblas
    ENDIF()
ENDIF()
IF (OPENBLAS_LIB)
    find_package( Threads REQUIRED )
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()
message( STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}" )


# Find libgfortran (static preferred)
IF(NOT WIN32)
    # Query gfortran to get the libgfortran path
    IF (CMAKE_Fortran_COMPILER)
        EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.a
                OUTPUT_VARIABLE _libgfortran_path
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
        IF (NOT EXISTS ${_libgfortran_path})
            EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.so
                    OUTPUT_VARIABLE _libgfortran_path
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
        ENDIF ()
    ENDIF()

    IF(EXISTS ${_libgfortran_path})
        get_filename_component(GFORTRAN_PATH ${_libgfortran_path} PATH)
        find_library(GFORTRAN_LIB gfortran PATHS ${GFORTRAN_PATH})
    ELSE()
        # if libgfortran wasn't found at this point, the installation is probably broken
        # Let's try to find the library nonetheless.
        FIND_LIBRARY(GFORTRAN_LIB gfortran)
    ENDIF()
    IF (NOT GFORTRAN_LIB)
        MESSAGE(STATUS "gfortran is required but it could not be found")
        SET(ABORT_CONFIG TRUE)
    ENDIF (NOT GFORTRAN_LIB)

    # we may need libquadmath.a , so add it if it is found
    find_library(QUADMATH_LIB quadmath PATHS ${GFORTRAN_PATH})
    IF (QUADMATH_LIB)
        list (APPEND GFORTRAN_LIB ${QUADMATH_LIB})
    ENDIF (QUADMATH_LIB)

    IF (APPLE)
        # also we need -lgcc_ext.10.5 on (most) MAC
        find_library(GCC_EXT_LIB gcc_ext.10.5 PATHS ${GFORTRAN_PATH})
        IF (GCC_EXT_LIB)
            list (APPEND GFORTRAN_LIB ${GCC_EXT_LIB})
        ELSE ()
            MESSAGE(STATUS "gcc_ext is required on MAC but could not be found")
            SET(ABORT_CONFIG TRUE)
        ENDIF ()
    ENDIF(APPLE)

ENDIF(NOT WIN32)


# ARPACK
find_library (ARPACK_LIB NAMES arpack PATHS external external/ARPACK)


IF (NOT ARPACK_LIB)
    MESSAGE(SEND_ERROR "libarpack is required but could not be found")
    SET(ABORT_CONFIG TRUE)
ENDIF()
message( STATUS "ARPACK_LIB found: ${ARPACK_LIB}" )

# SUPERLU
if (SUPERLU)
    find_path(SUPERLU_INC supermatrix.h "external/SuperLU/SRC" "/usr/local/include/superlu" )
    IF (NOT SUPERLU_INC)
        MESSAGE(STATUS "SuperLU (supermatrix.h) is required but could not be found")
        SET(ABORT_CONFIG TRUE)
    ENDIF()
    make_global_path(SUPERLU_INC "SUPERLU Include Path")

    find_library(SUPERLU_LIB superlu.a libsuperlu.a
            PATHS ${SUPERLU_INC} ${SUPERLU_INC}/../ external /usr/local/lib )
    IF (NOT SUPERLU_LIB)
        MESSAGE(STATUS "SuperLU Library is required but could not be found")
        SET(ABORT_CONFIG TRUE)
    ENDIF()
    make_global_path(SUPERLU_LIB "SUPERLU Library")
    include_directories( ${SUPERLU_INC} )
    message( STATUS "SUPERLU_LIB found: ${SUPERLU_LIB}" )
endif()

# Suitesparse
if (UMFPACK OR CHOLMOD)
    # Suitesparse CHOLMOD
    find_path(SUITESPARSE_DIR SuiteSparse_config/SuiteSparse_config.h "external/SuiteSparse")
    IF (NOT SUITESPARSE_DIR)
        # nothing in external, search globally for libraries and headers
        find_library(CHOLMOD_LIB libcholmod.a)
        find_library(COLAMD_LIB libcolamd.a)
        find_library(CCOLAMD_LIB libccolamd.a)
        find_library(METIS_LIB libmetis.a)
        IF (NOT METIS_LIB)
            find_library(METIS_LIB libmetis.dylib)
        ENDIF()
        IF (NOT METIS_LIB)
            find_library(METIS_LIB libmetis.so)
        ENDIF()
        find_library(CAMD_LIB libcamd.a)
        find_library(AMD_LIB libamd.a)

        # Suitesparse UMFPACK
        find_library(UMFPACK_LIB libumfpack.a)
        find_library(SUITESPARSE_LIB libsuitesparseconfig.a)

        FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
                PATHS /usr/local/include
                /usr/include
                /usr/include/suitesparse/
                ${CMAKE_SOURCE_DIR}/MacOS/Libs/cholmod
                PATH_SUFFIXES cholmod/ CHOLMOD/ )
        IF (NOT CHOLMOD_INCLUDE_DIR)
            MESSAGE(STATUS "SuiteSparse chomod.h is required but could not be found")
            SET(ABORT_CONFIG TRUE)
        ENDIF()

        include_directories(
                ${CHOLMOD_INCLUDE_DIR}
        )


    ELSE (NOT SUITESPARSE_DIR)
        # found suite sparse locally in external dir
        make_global_path(SUITESPARSE_DIR "SuiteSparse Directory")
        message( STATUS "SUITESPARSE_DIR found: ${SUITESPARSE_DIR}" )

        find_library(CHOLMOD_LIB libcholmod.a ${SUITESPARSE_DIR}/CHOLMOD/Lib  )
        find_library(COLAMD_LIB libcolamd.a ${SUITESPARSE_DIR}/COLAMD/Lib  )
        find_library(CCOLAMD_LIB libccolamd.a ${SUITESPARSE_DIR}/CCOLAMD/Lib  )
        find_library(METIS_LIB libmetis.a ${SUITESPARSE_DIR}/metis-4.0 )
        IF (NOT METIS_LIB)
            find_library(METIS_LIB libmetis.dylib)
        ENDIF()
        IF (NOT METIS_LIB)
            find_library(METIS_LIB libmetis.so)
        ENDIF()
        find_library(CAMD_LIB libcamd.a ${SUITESPARSE_DIR}/CAMD/Lib )
        find_library(AMD_LIB libamd.a ${SUITESPARSE_DIR}/AMD/Lib )

        # Suitesparse UMFPACK
        find_library(UMFPACK_LIB libumfpack.a ${SUITESPARSE_DIR}/UMFPACK/Lib  )
        find_library(SUITESPARSE_LIB libsuitesparseconfig.a ${SUITESPARSE_DIR}/SuiteSparse_config )

        include_directories(
                ${SUITESPARSE_DIR}/CHOLMOD/Include
                ${SUITESPARSE_DIR}/SuiteSparse_config
        )


    ENDIF()

    set(SSPARSE_LIBS
            ${UMFPACK_LIB}
            ${CHOLMOD_LIB}
            ${COLAMD_LIB}
            ${CCOLAMD_LIB}
            ${METIS_LIB}
            ${CAMD_LIB}
            ${AMD_LIB}
            ${SUITESPARSE_LIB}
            )

    # some linux systems need specifically librt (for SuiteSparse):
    if (NOT APPLE)
        set (SSPARSE_LIBS ${SSPARSE_LIBS} -lrt)
    endif()


endif()
add_executable(HMC HMC.cpp)
TARGET_LINK_LIBRARIES(HMC ${LAPACK_LIBRARIES})
TARGET_LINK_LIBRARIES(HMC ${ARPACK_LIB})
