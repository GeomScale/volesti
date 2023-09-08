set(LP_SOLVE_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetLPSolve)
  find_path(LP_SOLVE_DIR NAMES lpsolve.h PATHS ${LP_SOLVE_CMAKE_DIR}/../_deps/lpsolve-src)

  if (NOT LP_SOLVE_DIR)
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${LP_SOLVE_CMAKE_DIR}/../_deps")
      FetchContent_Declare(
          lpsolve
          URL https://webwerks.dl.sourceforge.net/project/lpsolve/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz
          URL_HASH MD5=a829a8d9c60ff81dc72ff52363703886
      )

      FetchContent_GetProperties(lpsolve)

      if(NOT lpsolve_POPULATED)
          message(STATUS "lp_solve library not found locally, downloading it.")
          FetchContent_Populate(lpsolve)
      endif()

      set(LP_SOLVE_DIR "${lpsolve_SOURCE_DIR}")
      message(STATUS "Using downloaded lp_solve at: ${LP_SOLVE_DIR}")

  else()

      message(STATUS "lp_solve library found: ${LP_SOLVE_DIR}")

  endif()

  include_directories(${LP_SOLVE_DIR})

endfunction()
