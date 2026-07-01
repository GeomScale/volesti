# Convex body tests

add_executable(convex_body_test convex_body_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME convex_body_test_basic COMMAND convex_body_test -tc=convex_body)
TARGET_LINK_LIBRARIES(convex_body_test coverage_config)

add_executable(order_polytope order_polytope.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME order_polytope_basics COMMAND order_polytope -tc=basics)
add_test(NAME order_polytope_line_intersect COMMAND order_polytope -tc=line_intersect)
add_test(NAME order_polytope_reflection COMMAND order_polytope -tc=reflection)
add_test(NAME order_polytope_vec_mult COMMAND order_polytope -tc=vec_mult)
TARGET_LINK_LIBRARIES(order_polytope lp_solve coverage_config)

add_subdirectory(correlation_matrices)
add_subdirectory(spectrahedra)
