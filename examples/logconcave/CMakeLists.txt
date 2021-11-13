# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.
# Licensed under GNU LGPL.3, see LICENCE file

option(USE_MKL "Use MKL library to build eigen" ON)
set(MKLROOT $ENV{HOME}/intel/mkl)

if (USE_MKL)
  find_library(BLAS NAME libblas.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  include_directories (BEFORE ${MKLROOT}/include)
  set(PROJECT_LIBS ${BLAS_LIBRARIES})
  add_definitions(-DEIGEN_USE_MKL_ALL)
endif(USE_MKL)

add_executable (simple_hmc simple_hmc.cpp)
add_executable (simple_ode simple_ode.cpp)
add_executable (simple_uld simple_uld.cpp)
add_executable (simple_sde simple_sde.cpp)
add_executable (high_dimensional_hmc high_dimensional_hmc.cpp)
add_executable (high_dimensional_hmc_rand_poly high_dimensional_hmc_rand_poly.cpp)
add_executable (exponential_exact_hmc exponential_exact_hmc.cpp)
add_executable (gaussian_exact_hmc gaussian_exact_hmc.cpp)

#add_executable (adaptive_step_size adaptive_step_size.cpp)

TARGET_LINK_LIBRARIES(high_dimensional_hmc_rand_poly ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
TARGET_LINK_LIBRARIES(high_dimensional_hmc ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
#TARGET_LINK_LIBRARIES(adaptive_step_size ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
TARGET_LINK_LIBRARIES(simple_ode ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
TARGET_LINK_LIBRARIES(simple_sde ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
TARGET_LINK_LIBRARIES(simple_uld ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
TARGET_LINK_LIBRARIES(simple_hmc ${LP_SOLVE} ${BLAS}  "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
