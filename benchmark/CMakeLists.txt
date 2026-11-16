CMAKE_MINIMUM_REQUIRED(VERSION 3.11)
project(VolestiBenchmark)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(VOLESTI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
endif(DISABLE_NLP_ORACLES)

include("${VOLESTI_ROOT}/external/cmake-files/Eigen.cmake")
GetEigen()

include("${VOLESTI_ROOT}/external/cmake-files/Boost.cmake")
GetBoost()
find_package(Boost REQUIRED COMPONENTS program_options)

include("${VOLESTI_ROOT}/external/cmake-files/HiGHS.cmake")
GetHighs()

include("${VOLESTI_ROOT}/external/cmake-files/QD.cmake")
GetQD()

include("${VOLESTI_ROOT}/external/cmake-files/Spectra.cmake")
GetSpectra()

set(HIGHS highs::highs)

if (NOT HIGHS)
  message(FATAL_ERROR "This program requires the HiGHS library, and will not be compiled.")
else ()
  message(STATUS "Library HiGHS found: ${HIGHS}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  include_directories (BEFORE ${VOLESTI_ROOT}/external)
  include_directories (BEFORE ${VOLESTI_ROOT}/include)
  include_directories (BEFORE ${VOLESTI_ROOT}/include/generators)
  include_directories (BEFORE ${VOLESTI_ROOT}/include/convex_bodies)
  include_directories (BEFORE ${VOLESTI_ROOT}/include/random_walks)
  include_directories (BEFORE ${VOLESTI_ROOT}/include/preprocess)

  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()

  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++17") 
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")         
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")

  add_executable (benchmark_run
    src/benchmark_main.cpp
    src/benchmark_run.cpp
    src/walk_parameters.cpp
    src/benchmark_utils.cpp
    src/walk_result.cpp
    src/walk_registry.cpp
  )

  TARGET_LINK_LIBRARIES(benchmark_run ${HIGHS} QD_LIB Boost::program_options)

endif()