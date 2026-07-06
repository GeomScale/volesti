# ODE solver tests

add_executable(ode_solvers_test ode_solvers_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME ode_solvers_test_first_order COMMAND ode_solvers_test -tc=first_order)
add_test(NAME ode_solvers_test_second_order COMMAND ode_solvers_test -tc=second_order)
TARGET_LINK_LIBRARIES(ode_solvers_test lp_solve ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${MKL_LINK} QD_LIB coverage_config)

add_executable(leapfrog_test leapfrog_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME leapfrog_test_no_boundaries COMMAND leapfrog_test -tc=generalized_leapfrog_no_boundaries)
add_test(NAME leapfrog_test_two_pairs COMMAND leapfrog_test -tc=generalized_leapfrog_two_pairs)
add_test(NAME leapfrog_test_constrained COMMAND leapfrog_test -tc=generalized_leapfrog_constrained)
add_test(NAME leapfrog_test_adaptive_flag COMMAND leapfrog_test -tc=generalized_leapfrog_adaptive_flag)
add_test(NAME leapfrog_test_get_set_state COMMAND leapfrog_test -tc=generalized_leapfrog_get_set_state)
add_test(NAME leapfrog_test_isotropic_quadratic COMMAND leapfrog_test -tc=generalized_leapfrog_isotropic_quadratic_functor)
add_test(NAME leapfrog_test_2d COMMAND leapfrog_test -tc=generalized_leapfrog_2d)
TARGET_LINK_LIBRARIES(leapfrog_test coverage_config)
