project( GSoC_2018_Rvolesti )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

add_executable (vol vol.cpp)

include_directories (BEFORE ../external/Eigen)
include_directories (BEFORE ../include/cartesian_geom)
include_directories (BEFORE ../include/convex_bodies)
include_directories (BEFORE ../include/LPsolve)
include_directories (BEFORE ../include)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/path/to/liblpsolve55.so)  # <-- edit the path to liblpsolve55.so
TARGET_LINK_LIBRARIES(vol liblpsolve55.so)
#add_library( liblpsolve55.so SHARED IMPORTED )
#TARGET_LINK_LIBRARIES(GLBall liblpaolve.so)
#add_library (liblpsolve55.so STATIC)  #build .so library

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
add_definitions(${CMAKE_CXX_FLAGS} "-lm")
add_definitions(${CMAKE_CXX_FLAGS} "-ldl")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl")




