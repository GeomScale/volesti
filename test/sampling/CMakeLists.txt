# MMCS sampling tests

add_executable(mmcs_test mmcs_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME test_mmcs COMMAND mmcs_test -tc=mmcs)
TARGET_LINK_LIBRARIES(mmcs_test lp_solve ${MKL_LINK} coverage_config)
