cmake_minimum_required(VERSION 3.17)

project(volesti LANGUAGES CXX)

option(USE_HIGHS "Use HiGHS for LP support when available" ON)

add_library(volesti INTERFACE)
target_include_directories(volesti INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(VOLESTI_LP_TARGET "")

if(USE_HIGHS)
  find_package(HiGHS QUIET)

  if(TARGET highs::highs)
    set(VOLESTI_LP_TARGET highs::highs)
  elseif(TARGET HiGHS::HiGHS)
    set(VOLESTI_LP_TARGET HiGHS::HiGHS)
  elseif(HiGHS_FOUND AND DEFINED HiGHS_LIBRARIES)
    set(VOLESTI_LP_TARGET "${HiGHS_LIBRARIES}")
    if(DEFINED HiGHS_INCLUDE_DIRS)
      target_include_directories(volesti INTERFACE ${HiGHS_INCLUDE_DIRS})
    endif()
  endif()
endif()

if(NOT VOLESTI_LP_TARGET)
  find_package(GLPK REQUIRED)

  if(TARGET GLPK::GLPK)
    set(VOLESTI_LP_TARGET GLPK::GLPK)
  elseif(DEFINED GLPK_LIBRARIES)
    set(VOLESTI_LP_TARGET "${GLPK_LIBRARIES}")
    if(DEFINED GLPK_INCLUDE_DIRS)
      target_include_directories(volesti INTERFACE ${GLPK_INCLUDE_DIRS})
    endif()
  else()
    message(FATAL_ERROR "GLPK was found but no usable CMake target or library variable is available")
  endif()
endif()

target_link_libraries(volesti INTERFACE ${VOLESTI_LP_TARGET})
