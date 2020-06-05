# VolEsti (volume computation and sampling library)
# Copyright (c) 20012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

enable_testing()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#if(NOT DEFINED BOOST)
#  message(FATAL_ERROR "This program requires the boost library, and will not be compiled.")
#else()
#  message(STATUS "BOOST is DEFINED")
#endif()

find_library(LP_SOLVE NAMES liblpsolve55.so liblpsolve55.dylib PATHS /usr/lib/lp_solve /opt/local/lib)

if (NOT LP_SOLVE)

  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")

else ()

  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../external/Eigen)
  include_directories (BEFORE ../external)
  include_directories (BEFORE ../external/minimum_ellipsoid)
  #include_directories (BEFORE ../include/cartesian_geom)
  #include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../external/LPsolve_src/run_headers)
  include_directories (BEFORE ../external/boost)
  #include_directories (BEFORE BOOST)
  include_directories (BEFORE ../include/generators)
  include_directories (BEFORE ../include/volume)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../include/convex_bodies/spectrahedra)
  include_directories (BEFORE ../include/annealing)
  include_directories (BEFORE ../include/samplers)
  include_directories (BEFORE ../include/lp_oracles)
  include_directories (BEFORE ../include/misc)

  #for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0") 
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()
  
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-g")  # enable debuger
  #add_definitions(${CMAKE_CXX_FLAGS} "-Wint-in-bool-context")
  #add_definitions(${CMAKE_CXX_FLAGS} "-Wall")

  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )
 
  add_executable (new_volume_example new_volume_example.cpp)
  add_executable (benchmarks_sob benchmarks_sob.cpp)
  add_executable (benchmarks_cg benchmarks_cg.cpp)
  add_executable (benchmarks_cb benchmarks_cb.cpp)

  add_library(test_main OBJECT test_main.cpp)


  add_executable(test_sdpa_format test_sdpa_format.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME test_sdpa_format COMMAND test_sdpa_format -tc=sdpa_format_parser)

  add_executable (volume_sob_hpolytope volume_sob_hpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_sob_hpolytope_cube COMMAND volume_sob_hpolytope -tc=cube)
  add_test(NAME volume_sob_hpolytope_cross COMMAND volume_sob_hpolytope -tc=cross)
  add_test(NAME volume_sob_hpolytope_birkhoff COMMAND volume_sob_hpolytope -tc=birk)
  add_test(NAME volume_sob_hpolytope_prod_simplex COMMAND volume_sob_hpolytope -tc=prod_simplex)
  add_test(NAME volume_sob_hpolytope_simplex COMMAND volume_sob_hpolytope -tc=simplex)
  add_test(NAME volume_sob_hpolytope_skinny_cube COMMAND volume_sob_hpolytope -tc=skinny_cube)

  add_executable (volume_sob_vpolytope volume_sob_vpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_sob_vpolytope_cube COMMAND volume_sob_vpolytope -tc=cube)
  add_test(NAME volume_sob_vpolytope_cross COMMAND volume_sob_vpolytope -tc=cross)
  add_test(NAME volume_sob_vpolytope_simplex COMMAND volume_sob_vpolytope -tc=simplex)

  add_executable (volume_cg_hpolytope volume_cg_hpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cg_hpolytope_cube COMMAND volume_cg_hpolytope -tc=cube)
  add_test(NAME volume_cg_hpolytope_cross COMMAND volume_cg_hpolytope -tc=cross)
  add_test(NAME volume_cg_hpolytope_birkhoff COMMAND volume_cg_hpolytope -tc=birk)
  add_test(NAME volume_cg_hpolytope_prod_simplex COMMAND volume_cg_hpolytope -tc=prod_simplex)
  add_test(NAME volume_cg_hpolytope_simplex COMMAND volume_cg_hpolytope -tc=simplex)
  add_test(NAME volume_cg_hpolytope_skinny_cube COMMAND volume_cg_hpolytope -tc=skinny_cube)

  add_executable (volume_cg_vpolytope volume_cg_vpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cg_vpolytope_cube COMMAND volume_cg_vpolytope -tc=cube)
  add_test(NAME volume_cg_vpolytope_cross COMMAND volume_cg_vpolytope -tc=cross)
  add_test(NAME volume_cg_vpolytope_simplex COMMAND volume_cg_vpolytope -tc=simplex)

  add_executable (volume_cb_hpolytope volume_cb_hpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_hpolytope_cube COMMAND volume_cb_hpolytope -tc=cube)
  add_test(NAME volume_cb_hpolytope_cross COMMAND volume_cb_hpolytope -tc=cross)
  add_test(NAME volume_cb_hpolytope_birkhoff COMMAND volume_cb_hpolytope -tc=birk)
  add_test(NAME volume_cb_hpolytope_prod_simplex COMMAND volume_cb_hpolytope -tc=prod_simplex)
  add_test(NAME volume_cb_hpolytope_simplex COMMAND volume_cb_hpolytope -tc=simplex)
  add_test(NAME volume_cb_hpolytope_skinny_cube COMMAND volume_cb_hpolytope -tc=skinny_cube)

  add_executable (volume_cb_vpolytope volume_cb_vpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_vpolytope_cube COMMAND volume_cb_vpolytope -tc=cube)
  add_test(NAME volume_cb_vpolytope_cross COMMAND volume_cb_vpolytope -tc=cross)
  add_test(NAME volume_cb_vpolytope_simplex COMMAND volume_cb_vpolytope -tc=simplex)

  add_executable (volume_cb_zonotopes volume_cb_zonotopes.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_zonotopes_uniform_zonotopes
           COMMAND volume_cb_zonotopes -tc=uniform_zonotopes)

  add_executable (volume_cb_vpoly_intersection_vpoly volume_cb_vpoly_intersection_vpoly.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_vpoly_intersection_vpoly_random_vpoly_sphere
           COMMAND volume_cb_vpoly_intersection_vpoly -tc=random_vpoly_sphere)

  add_executable (new_rounding_test new_rounding_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME new_rounding_test_round_skinny_cube
           COMMAND new_rounding_test -tc=round_skinny_cube)

  TARGET_LINK_LIBRARIES(new_volume_example ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_sob_hpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_sob_vpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cg_hpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cg_vpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_hpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_vpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_zonotopes ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_vpoly_intersection_vpoly ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(new_rounding_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(benchmarks_sob ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(benchmarks_cg ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(benchmarks_cb ${LP_SOLVE})

endif()


