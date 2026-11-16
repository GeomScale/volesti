
cmake_minimum_required(VERSION 3.11)
project(VolEstiGeneralSampling)

option(DISABLE_NLP_ORACLES "Disable non-linear oracles" ON)
option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3 -march=native -DTIME_KEEPING -DBOOST_NO_AUTO_PTR -DEIGEN_NO_DEBUG)

if(DISABLE_NLP_ORACLES)
    add_definitions(-DDISABLE_NLP_ORACLES)
endif()

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/HiGHS.cmake")
GetHighs()

include("../../external/cmake-files/QD.cmake")
GetQD()

set(HIGHS highs::highs)

if(NOT TARGET highs::highs)
    message(FATAL_ERROR "HiGHS target was not found!")
else()
    message(STATUS "Library HiGHS found: ${HIGHS}")
endif()

include_directories(BEFORE ../../include)
include_directories(BEFORE ../../include/convex_bodies)
include_directories(BEFORE ../../include/generators)
include_directories(BEFORE ../../external)
include_directories(BEFORE ../../external/minimum_ellipsoid)

find_package(Threads REQUIRED)

add_executable(general_sampling_time_limit
    general_sampling_time_limit.cpp
)

target_link_libraries(general_sampling_time_limit
    PRIVATE
    ${HIGHS}
    qd
    Threads::Threads
)
