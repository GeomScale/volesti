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

find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)

  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")

else ()

  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  #include_directories (BEFORE ../external/Eigen)
  include_directories (BEFORE ../external)
  include_directories (BEFORE ../external/minimum_ellipsoid)
  #include_directories (BEFORE ../include/cartesian_geom)
  #include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../external/LPsolve_src)
  include_directories (BEFORE ../external/boost)
  include_directories (BEFORE ../external/number_types)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../include/samplers)


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )

  add_executable (vol vol.cpp)

  add_library(test_main OBJECT test_main.cpp)

  add_executable (volume_test volume_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (cheb_test chebychev_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (rounding_test rounding_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (volumeCV_test volumeCV_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (long_double_CV long_double_testCV.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (float_VE float_testVE.cpp $<TARGET_OBJECTS:test_main>)
  
  add_test(NAME volume_cube COMMAND volume_test -tc=cube)
  add_test(NAME volume_cross COMMAND volume_test -tc=cross)
  add_test(NAME volume_birkhoff COMMAND volume_test -tc=birkhoff)
  add_test(NAME volume_prod_simplex COMMAND volume_test -tc=prod_simplex)
  add_test(NAME volume_skinny_cube COMMAND volume_test -tc=skinny_cube)

  add_test(NAME volumeCV_cube COMMAND volumeCV_test -tc=cube)
  add_test(NAME volumeCV_cross COMMAND volumeCV_test -tc=cross)
  add_test(NAME volumeCV_birkhoff COMMAND volumeCV_test -tc=birkhoff)
  add_test(NAME volumeCV_prod_simplex COMMAND volumeCV_test -tc=prod_simplex)

  add_test(NAME volume_float_cube COMMAND float_VE -tc=cube)
  add_test(NAME volume_float_cross COMMAND float_VE -tc=cross)
  add_test(NAME volume_float_birkhoff COMMAND float_VE -tc=birkhoff)
  add_test(NAME volume_float_prod_simplex COMMAND float_VE -tc=prod_simplex)
  add_test(NAME volume_float_skinny_cube COMMAND float_VE -tc=skinny_cube)

  add_test(NAME volume_long_double_cube COMMAND long_double_CV -tc=cube)
  add_test(NAME volumeCV_long_double_cross COMMAND long_double_CV -tc=cross)
  add_test(NAME volumeCV_long_double_birkhoff COMMAND long_double_CV -tc=birkhoff)
  add_test(NAME volumeCV_long_double_prod_simplex COMMAND long_double_CV -tc=prod_simplex)

  add_test(NAME cheb_cube COMMAND cheb_test -tc=cheb_cube)
  add_test(NAME cheb_cross COMMAND cheb_test -tc=cheb_cross)
  add_test(NAME cheb_birkhoff COMMAND cheb_test -tc=cheb_birkhoff)
  add_test(NAME cheb_prod_simplex COMMAND cheb_test -tc=cheb_prod_simplex)
  add_test(NAME cheb_simplex COMMAND cheb_test -tc=cheb_simplex)
  add_test(NAME cheb_skinny_cube COMMAND cheb_test -tc=cheb_skinny_cube)

  add_test(NAME round_skinny_cube COMMAND rounding_test -tc=round_skinny_cube)
  add_test(NAME round_rot_skinny_cube COMMAND rounding_test -tc=round_rot_skinny_cube)


  TARGET_LINK_LIBRARIES(vol ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(cheb_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(rounding_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volumeCV_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(float_VE ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(long_double_CV ${LP_SOLVE})

endif()



