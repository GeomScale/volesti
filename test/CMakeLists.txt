# VolEsti (volume computation and sampling library)
# Copyright (c) 20012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
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

find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

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

  #add_executable (vol vol.cpp)
  #add_executable (volume volume_example.cpp)
  #add_executable (generate generator.cpp)
 
  add_executable (new_volume_example new_volume_example.cpp)
  add_executable (benchmarks_sob benchmarks_sob.cpp)
  add_executable (benchmarks_cg benchmarks_cg.cpp)
  add_executable (benchmarks_cb benchmarks_cb.cpp)

  add_library(test_main OBJECT test_main.cpp)

  #add_executable (volume_test volume_test.cpp $<TARGET_OBJECTS:test_main>)

  #add_executable (cheb_test chebychev_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (rounding_test rounding_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (volumeCG_test volumeCG_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (VpolyCG_test VpolyCG_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (VpolyVol_test VpolyVol_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (ZonotopeVol_test ZonotopeVol_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (cool_bodies_bill_test cooling_bodies_bill_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (ZonotopeVolCG_test ZonotopeVolCG_test.cpp $<TARGET_OBJECTS:test_main>)

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

  #add_test(NAME volume_cube COMMAND volume_test -tc=cube)
  #add_test(NAME volume_cross COMMAND volume_test -tc=cross)
  #add_test(NAME volume_birkhoff COMMAND volume_test -tc=birk)
  #add_test(NAME volume_prod_simplex COMMAND volume_test -tc=prod_simplex)
  #add_test(NAME volume_simplex COMMAND volume_test -tc=simplex)
  #add_test(NAME volume_skinny_cube COMMAND volume_test -tc=skinny_cube)

  #add_test(NAME volumeCG_cube COMMAND volumeCG_test -tc=cube)
  #add_test(NAME volumeCG_cross COMMAND volumeCG_test -tc=cross)
  #add_test(NAME volumeCG_birkhoff COMMAND volumeCG_test -tc=birk)
  #add_test(NAME volumeCG_prod_simplex COMMAND volumeCG_test -tc=prod_simplex)
  #add_test(NAME volumeCG_simplex COMMAND volumeCG_test -tc=simplex)

  #add_test(NAME VpolyCG_cube COMMAND VpolyCG_test -tc=cube)
  #add_test(NAME VpolyCG_cross COMMAND VpolyCG_test -tc=cross)
  #add_test(NAME VpolyCG_simplex COMMAND VpolyCG_test -tc=simplex)
  
  #add_test(NAME VpolyVol_cube COMMAND VpolyVol_test -tc=cube)
  #add_test(NAME VpolyVol_cross COMMAND VpolyVol_test -tc=cross)
  #add_test(NAME VpolyVol_simplex COMMAND VpolyVol_test -tc=simplex)

  #add_test(NAME ZonotopeVol4 COMMAND ZonotopeVol_test -tc=4_dimensional)

  #add_test(NAME ZonotopeVolCG4 COMMAND ZonotopeVolCG_test -tc=4_dimensional)
  #add_test(NAME ZonotopeVolCG5 COMMAND ZonotopeVolCG_test -tc=5_dimensional)

  #add_test(NAME cheb_cube COMMAND cheb_test -tc=cheb_cube)
  #add_test(NAME cheb_cross COMMAND cheb_test -tc=cheb_cross)
  #add_test(NAME cheb_birkhoff COMMAND cheb_test -tc=cheb_birkhoff)
  #add_test(NAME cheb_prod_simplex COMMAND cheb_test -tc=cheb_prod_simplex)
  #add_test(NAME cheb_simplex COMMAND cheb_test -tc=cheb_simplex)
  #add_test(NAME cheb_skinny_cube COMMAND cheb_test -tc=cheb_skinny_cube)

  #add_test(NAME cool_bodies_cube COMMAND cool_bodies_bill_test -tc=cube)
  #add_test(NAME cool_bodies_cross COMMAND cool_bodies_bill_test -tc=cross)
  #add_test(NAME cool_bodies_birkhoff COMMAND cool_bodies_bill_test -tc=birk)
  #add_test(NAME cool_bodies_prod_simplex COMMAND cool_bodies_bill_test -tc=prod_simplex)
  #add_test(NAME cool_bodies_simplex COMMAND cool_bodies_bill_test -tc=simplex)
  #add_test(NAME cool_bodies_skinny_cube COMMAND cool_bodies_bill_test -tc=skinny_cube)

  #add_test(NAME round_skinny_cube COMMAND rounding_test -tc=round_skinny_cube)
  #add_test(NAME round_rot_skinny_cube COMMAND rounding_test -tc=round_rot_skinny_cube)

  #TARGET_LINK_LIBRARIES(vol ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volume ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(generate ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volume_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(cheb_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(rounding_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volumeCG_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(VpolyCG_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(VpolyVol_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(ZonotopeVol_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(cool_bodies_bill_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(ZonotopeVolCG_test ${LP_SOLVE})
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


