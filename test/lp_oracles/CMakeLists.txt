# LP oracle and boundary oracle tests

add_executable(boundary_oracles_test boundary_oracles_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME boundary_oracles_test_h_poly_oracles COMMAND boundary_oracles_test -tc=h_poly_oracles)
TARGET_LINK_LIBRARIES(boundary_oracles_test lp_solve ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${MKL_LINK} coverage_config)
