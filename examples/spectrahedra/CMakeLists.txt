# VolEsti (volume computation and sampling library)
# Copyright (c) 20012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
# Licensed under GNU LGPL.3, see LICENCE file


add_executable (read_write_sdpa_file read_write_sdpa_file.cpp)


# Find LAPACK and BLAS
# OPENBLAS or ( ( SystemOpenblas or BLAS) and LAPACK)
## prefer local openblas
find_library(OPENBLAS_LIB openblas PATHS external NO_DEFAULT_PATH)
IF (OPENBLAS_LIB)
    set(LAPACK_LIBRARIES ${OPENBLAS_LIB}) #local openblas has lapack build in

    find_package( Threads REQUIRED )
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    message( STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}" )


    # ARPACK
    find_library(ARPACK_LIB arpack PATHS external NO_DEFAULT_PATH)
    IF (ARPACK_LIB)
        message( STATUS "ARPACK_LIB found: ${ARPACK_LIB}" )

        # Query gfortran to get the libgfortran path
        FIND_LIBRARY(GFORTRAN_LIB NAMES libgfortran.so PATHS /usr/lib/gcc/x86_64-linux-gnu/8/)

        IF (GFORTRAN_LIB)
            message( STATUS "GFORTRAN_LIB found: ${GFORTRAN_LIB}" )


            add_executable (boltzmann_hmc_walk boltzmann_hmc_walk.cpp)
            TARGET_LINK_LIBRARIES(boltzmann_hmc_walk ${GFORTRAN_LIB} ${ARPACK_LIB} ${LAPACK_LIBRARIES})

            add_executable (solve_sdp solve_sdp.cpp)
            TARGET_LINK_LIBRARIES(solve_sdp ${GFORTRAN_LIB} ${ARPACK_LIB} ${LAPACK_LIBRARIES})

        ELSE()
            MESSAGE(STATUS "gfortran is required but it could not be found")
        ENDIF ()


    ELSE()
        message(FATAL_ERROR "This program requires the arpack library, and will not be compiled.")
    ENDIF()

ELSE()
    message(FATAL_ERROR "This program requires the openblas library, and will not be compiled.")
ENDIF()






