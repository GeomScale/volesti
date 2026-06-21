cmake_minimum_required(VERSION 3.10)
project(VolestiBenchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3)
add_compile_definitions(EIGEN_NO_DEBUG)
add_compile_definitions(BOOST_NO_AUTO_PTR)
add_compile_definitions(DISABLE_NLP_ORACLES) 

find_package(Boost REQUIRED COMPONENTS program_options)
find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve /usr/lib /usr/local/lib)
find_library(QD_LIB NAMES qd PATHS /usr/lib /usr/local/lib)

if (NOT LP_SOLVE)
    message(WARNING "lp_solve not found. ComputeInnerBall() might fail to link.")
else()
    message(STATUS "Library lp_solve found: ${LP_SOLVE}")
endif()

if (NOT QD_LIB)
    message(WARNING "qd library not found. CRHMC / dd_real math will fail to link.")
else()
    message(STATUS "Library qd found: ${QD_LIB}")
endif()

include_directories(/usr/include/eigen3)
include_directories(/usr/include/lpsolve)
include_directories(/usr/include/lp_solve)
include_directories(include)
include_directories(../include) 
include_directories(../external)

include_directories(../include) 
include_directories(../include/generators)
include_directories(../include/convex_bodies)
include_directories(../include/random_walks)
include_directories(../include/preprocess)

set(SOURCES
    src/benchmark_main.cpp
    src/walk_parameters.cpp
    src/benchmark_utils.cpp
    src/walk_result.cpp
    src/walk_registry.cpp
)

add_executable(benchmark_run ${SOURCES})

if (LP_SOLVE)
    target_link_libraries(benchmark_run Boost::program_options ${LP_SOLVE} ${QD_LIB} m)
else()
    target_link_libraries(benchmark_run Boost::program_options ${QD_LIB} m)
endif()