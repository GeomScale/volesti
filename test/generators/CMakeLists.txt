# Polytope generator tests

add_executable(generators_test generators_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME test_known_polytope_generators COMMAND generators_test -tc=known_polytope_generators)
TARGET_LINK_LIBRARIES(generators_test lp_solve ${MKL_LINK} coverage_config)

add_executable(zonotope_gen_test zonotope_gen_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME zonotope_gen_test COMMAND zonotope_gen_test -tc=zonotope_generators)
TARGET_LINK_LIBRARIES(zonotope_gen_test lp_solve coverage_config)
