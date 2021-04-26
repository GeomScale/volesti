function(GetLPsolver)
    find_path(LPSOLVE_DIR NAMES LPsolve_src PATHS ../../external)

  if (NOT LPSOLVE_DIR) 
    include(FetchContent)
    FetchContent_Declare(
        lpsolve
        GIT_REPOSITORY https://github.com/gaborcsardi/lpSolve.git 
        GIT_TAG 5.6.13.3
    )

    FetchContent_GetProperties(lpSolve)
    
    if(NOT lpsolve_POPULATED)
      message(STATUS "lpSolve not found locally, downloading it.")
      FetchContent_Populate(lpsolve)
    endif()

    set(LPSOLVE_DIR "${lpsolve_SOURCE_DIR}/src")
    message(STATUS "Using downloaded lpSolve library at: ${LPSOLVE_DIR}")

  else ()

    message(STATUS "Eigen Library found: ${LPSOLVE_DIR}")

endif()

  include_directories(${LPSOLVE_DIR})

endfunction()
