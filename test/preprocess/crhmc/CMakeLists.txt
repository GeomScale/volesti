# CRHMC preprocessing tests

set(ADDITIONAL_FLAGS "-march=native -DSIMD_LEN=0 -DTIME_KEEPING")

add_executable(crhmc_polytope_preparation_test crhmc_polytope_preparation_test.cpp $<TARGET_OBJECTS:test_main>)
# test_preparation_crhmc loads e_coli metabolic network and takes 100+ min — run manually only
# add_test(NAME crhmc_polytope_test_preparation COMMAND crhmc_polytope_preparation_test -tc=test_preparation_crhmc)
add_test(NAME crhmc_test_fixed_vars COMMAND crhmc_polytope_preparation_test -tc=test_fixed_vars_crhmc)
add_test(NAME crhmc_test_dep_vars COMMAND crhmc_polytope_preparation_test -tc=test_dep_vars_crhmc)
add_test(NAME crhmc_test_center_computation COMMAND crhmc_polytope_preparation_test -tc=test_center_computation)
set_target_properties(crhmc_polytope_preparation_test PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS})
TARGET_LINK_LIBRARIES(crhmc_polytope_preparation_test ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${MKL_LINK} QD_LIB coverage_config)
