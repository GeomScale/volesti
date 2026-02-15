cmake_minimum_required(VERSION 3.11)
project(VolEstiGeneralSampling)

# --------------------------------------------------------
# 1. Options & Flags
# --------------------------------------------------------
option(DISABLE_NLP_ORACLES "Disable non-linear oracles" ON)
option(BUILTIN_EIGEN "Use eigen from ../../external" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization & Standard Flags
add_compile_options(-O3 -march=native -DTIME_KEEPING -DBOOST_NO_AUTO_PTR -DEIGEN_NO_DEBUG)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
endif()

# --------------------------------------------------------
# 2. Dependencies (Using Relative Paths)
# --------------------------------------------------------
# We use the helper scripts located in the main library
include("../../external/cmake-files/Eigen.cmake")
GetEigen()

include("../../external/cmake-files/Boost.cmake")
GetBoost()

include("../../external/cmake-files/LPSolve.cmake")
GetLPSolve()

include("../../external/cmake-files/QD.cmake")
GetQD() # This sets the variable ${QD_LIB}

# --------------------------------------------------------
# 3. Find Libraries Manually
# --------------------------------------------------------
# Find lp_solve (standard paths + typical linux paths)
find_library(LP_SOLVE NAMES liblpsolve55.so lpsolve55 PATHS 
    /usr/lib/lp_solve 
    /usr/lib/x86_64-linux-gnu 
    /usr/local/lib
)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "lp_solve library not found! Please install it.")
else ()
  message(STATUS "Library lp_solve found: ${LP_SOLVE}")
endif()

# --------------------------------------------------------
# 4. Include Directories 
# --------------------------------------------------------
# This gives you access to "volesti/include/..."
include_directories(BEFORE ../../include)
include_directories(BEFORE ../../external)
include_directories(BEFORE ../../external/minimum_ellipsoid)

# --------------------------------------------------------
# 5. Targets
# --------------------------------------------------------
add_executable(general_sampling_time_limit general_sampling_time_limit.cpp)
target_link_libraries(general_sampling_time_limit ${LP_SOLVE} qd pthread)
