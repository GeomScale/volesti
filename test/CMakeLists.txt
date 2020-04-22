# VolEsti (volume computation and sampling library)
# Copyright (c) 20012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

#enable_testing()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

#if (NOT LP_SOLVE)

#  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")

#else ()

  #message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../external/Eigen)
  include_directories (BEFORE ../external/Spectra)
  include_directories (BEFORE ../external)
  include_directories (BEFORE ../external/minimum_ellipsoid)
  include_directories (BEFORE ../external/unsupported)
  #include_directories (BEFORE ../include/cartesian_geom)
  #include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../external/LPsolve_src/run_headers)
  include_directories (BEFORE ../external/boost)
  include_directories (BEFORE ../include/generators)
  include_directories (BEFORE ../include/volume)
  include_directories (BEFORE ../include/optimization)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../include/annealing)
  include_directories (BEFORE ../include/samplers)


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-march=native")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-msse2")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )

  add_executable (vol vol.cpp)
  add_executable (generate generate_sdp.cpp)
  #add_executable (generate_lp generate_lp.cpp)
  #add_executable (generate_sdp generate_sdp.cpp)

  
  


  #add_test(NAME round_skinny_cube COMMAND rounding_test -tc=round_skinny_cube)
  #add_test(NAME round_rot_skinny_cube COMMAND rounding_test -tc=round_rot_skinny_cube)

  #TARGET_LINK_LIBRARIES(vol ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(opt_cutting_plane ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(opt_cutting_plane_sdp ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volume ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(generate ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volume_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(cheb_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(rounding_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volumeCV_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(VpolyCV_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(VpolyVol_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(ZonotopeVol_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(ZonotopeVolCV_test ${LP_SOLVE})

#endif()


