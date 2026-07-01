# Monte Carlo integration tests

add_executable(simple_mc_integration simple_mc_integration.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME simple_mc_integration_over_limits COMMAND simple_mc_integration -tc=rectangle)
add_test(NAME simple_mc_integration_over_cubes COMMAND simple_mc_integration -tc=cube)
add_test(NAME simple_mc_integration_over_simplices COMMAND simple_mc_integration -tc=simplex)
add_test(NAME simple_mc_integration_over_product_simplices COMMAND simple_mc_integration -tc=prod_simplex)
add_test(NAME simple_mc_integration_over_cross_polytopes COMMAND simple_mc_integration -tc=cross)
add_test(NAME simple_mc_integration_over_birkhoff_polytopes COMMAND simple_mc_integration -tc=birkhoff)
TARGET_LINK_LIBRARIES(simple_mc_integration lp_solve ${MKL_LINK} coverage_config)
