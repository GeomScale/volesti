cmake_minimum_required(VERSION 3.12)
project(LP_IPM)

set(CMAKE_CXX_STANDARD 20)

add_executable(LP_IPM main.cpp InteriorPointMethod.cpp NewtonsMethod.cpp)

target_link_libraries(LP_IPM armadillo)
