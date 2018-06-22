project( VolEsti )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

find_library(LP_SOLVE NAMES liblpsolve55.so)

if (LP_SOLVE-NOTFOUND)

  message(STATUS "This program requires the lp_solve library, and will not be compiled.")

else ()

  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../external/Eigen)
  include_directories (BEFORE ../include/cartesian_geom)
  include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../include/LPsolve)
  include_directories (BEFORE ../include)

  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-lm")
  add_definitions(${CMAKE_CXX_FLAGS} "-ldl")

  add_executable (vol vol.cpp)

  TARGET_LINK_LIBRARIES(vol ${LP_SOLVE})

endif()



