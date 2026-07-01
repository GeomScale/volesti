# CRHMC sampling tests

set(ADDITIONAL_FLAGS "-march=native -DSIMD_LEN=0 -DTIME_KEEPING")

add_executable(crhmc_sampling_test crhmc_sampling_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME crhmc_sampling_test_crhmc COMMAND crhmc_sampling_test -tc=crhmc)
add_test(NAME crhmc_test_polytope_sampling COMMAND crhmc_sampling_test -tc=test_polytope_sampling_crhmc)
add_test(NAME crhmc_test_sparse_sampling COMMAND crhmc_sampling_test -tc=test_sampling_sparse_problem)
set_target_properties(crhmc_sampling_test PROPERTIES COMPILE_FLAGS ${ADDITIONAL_FLAGS})
TARGET_LINK_LIBRARIES(crhmc_sampling_test lp_solve ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${MKL_LINK} QD_LIB coverage_config)
