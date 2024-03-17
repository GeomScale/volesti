project( VolEsti-cpp-example )

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

add_definitions(-DDISABLE_NLP_ORACLES)

include("../../external/cmake-files/Eigen.cmake")
GetEigen()
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  add_compile_options(-D "EIGEN_NO_DEBUG")
else ()
  add_compile_definitions("EIGEN_NO_DEBUG")
endif ()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

include_directories (BEFORE ../../external)
include_directories (BEFORE ../../include)

add_executable (volume_example volume_example.cpp)
target_link_libraries(volume_example PUBLIC lp_solve)