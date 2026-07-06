# Root finder tests

add_executable(root_finders_test root_finders_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME root_finders_test_root_finders COMMAND root_finders_test -tc=root_finders)
TARGET_LINK_LIBRARIES(root_finders_test ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${MKL_LINK} coverage_config)
