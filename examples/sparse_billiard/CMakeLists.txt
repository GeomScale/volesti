cmake_minimum_required(VERSION 3.10)
project(SparseBilliardWalkBenchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
endif()

include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS random)

include_directories(${CMAKE_SOURCE_DIR}/../../include)
include_directories(${CMAKE_SOURCE_DIR}/../../external)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(sparse_billiard_benchmark sparse_billiard_walk_benchmark.cpp)

target_link_libraries(sparse_billiard_benchmark 
    lp_solve
    ${Boost_LIBRARIES}
    Eigen3::Eigen
)

target_compile_options(sparse_billiard_benchmark PRIVATE 
    -O3 
    -DEIGEN_NO_DEBUG
    -march=native
    -funroll-loops
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sparse_billiard_benchmark OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(sparse_billiard_benchmark PRIVATE DISABLE_NLP_ORACLES) 