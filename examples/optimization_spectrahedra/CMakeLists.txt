# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
# Licensed under GNU LGPL.3, see LICENCE file

cmake_minimum_required(VERSION 3.16)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


add_executable (read_write_sdpa_file read_write_sdpa_file.cpp)

include("Eigen.cmake")
GetEigen()

include("Boost.cmake")
GetBoost()

include_directories (BEFORE ../../external)
include_directories (BEFORE ../../include/generators)
include_directories (BEFORE ../../include/volume)
include_directories (BEFORE ../../include)
include_directories (BEFORE ../../include/convex_bodies)
include_directories (BEFORE ../../include/convex_bodies/spectrahedra)
include_directories (BEFORE ../../include/annealing)
include_directories (BEFORE ../../include/samplers)
include_directories (BEFORE ../../include/misc)

# Find LAPACK and BLAS
# OPENBLAS or ( ( SystemOpenblas or BLAS) and LAPACK)
## prefer local openblas
find_library(OPENBLAS_LIB openblas PATHS external /usr/lib/gcc/x86_64-linux-gnu/)
IF (OPENBLAS_LIB)
    set(LAPACK_LIBRARIES ${OPENBLAS_LIB}) #local openblas has lapack build in

    find_package( Threads REQUIRED )
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    message( STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}" )


    # ARPACK
    find_library(ARPACK_LIB arpack PATHS external /usr/lib/gcc/x86_64-linux-gnu/)
    IF (ARPACK_LIB)
        message( STATUS "ARPACK_LIB found: ${ARPACK_LIB}" )

        # Query gfortran to get the libgfortran path
        FIND_LIBRARY(GFORTRAN_LIB NAMES libgfortran.so PATHS /usr/lib/gcc/x86_64-linux-gnu/9/)

        # Find libgfortran (static preferred)
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



        IF (GFORTRAN_LIB)
            message( STATUS "GFORTRAN_LIB found: ${GFORTRAN_LIB}" )

            # These are not compiling, see issue https://github.com/GeomScale/volesti/issues/264
            #add_executable (boltzmann_hmc_walk boltzmann_hmc_walk.cpp)
            #TARGET_LINK_LIBRARIES(boltzmann_hmc_walk ${ARPACK_LIB} ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})

            #add_executable (solve_sdp solve_sdp.cpp)
            #TARGET_LINK_LIBRARIES(solve_sdp ${GFORTRAN_LIB} ${LAPACK_LIBRARIES} ${ARPACK_LIB})

        ELSE()
            MESSAGE(STATUS "gfortran is required but it could not be found")
        ENDIF ()


    ELSE()
        message(FATAL_ERROR "This program requires the arpack library, and will not be compiled.")
    ENDIF()

ELSE()
    message(FATAL_ERROR "This program requires the openblas library, and will not be compiled.")
ENDIF()
