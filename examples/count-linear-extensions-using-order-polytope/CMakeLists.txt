# VolEsti (volume computation and sampling library)

add_executable (volesti_lecount_order_polytope volesti_lecount_order_polytope.cpp)
TARGET_LINK_LIBRARIES(volesti_lecount_order_polytope ${LP_SOLVE})
