cmake_minimum_required(VERSION 3.10)
project(VolestiBenchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-O3)

add_compile_definitions(EIGEN_NO_DEBUG)
add_compile_definitions(BOOST_NO_AUTO_PTR)
add_compile_definitions(DISABLE_NLP_ORACLES)

find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_library(LP_SOLVE NAMES liblpsolve55.so
             PATHS /usr/lib/lp_solve /usr/local/lib)
find_path(LP_SOLVE_INCLUDE_DIR NAMES lp_lib.h
          PATH_SUFFIXES lpsolve lp_solve)

find_library(QD_LIB NAMES qd)

if(NOT LP_SOLVE)
    message(WARNING "lp_solve not found. ComputeInnerBall() might fail to link.")
else()
    message(STATUS "Library lp_solve found: ${LP_SOLVE}")
endif()

if(NOT LP_SOLVE_INCLUDE_DIR)
    message(WARNING "lp_solve headers not found (lp_lib.h). Adjust include paths if needed.")
else()
    message(STATUS "lp_solve headers found: ${LP_SOLVE_INCLUDE_DIR}")
endif()

if(NOT QD_LIB)
    message(WARNING "qd library not found. CRHMC / dd_real math will fail to link.")
else()
    message(STATUS "Library qd found: ${QD_LIB}")
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../external
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/generators
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/convex_bodies
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/random_walks
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/preprocess
)

if(LP_SOLVE_INCLUDE_DIR)
    include_directories(${LP_SOLVE_INCLUDE_DIR})
endif()

set(SOURCES
    src/benchmark_main.cpp
    src/benchmark_run.cpp
    src/walk_parameters.cpp
    src/benchmark_utils.cpp
    src/walk_result.cpp
    src/walk_registry.cpp
)

add_executable(benchmark_run ${SOURCES})

target_link_libraries(benchmark_run
    Boost::program_options
    Eigen3::Eigen
    m
)

if(LP_SOLVE)
    target_link_libraries(benchmark_run ${LP_SOLVE})
endif()

if(QD_LIB)
    target_link_libraries(benchmark_run ${QD_LIB})
endif()