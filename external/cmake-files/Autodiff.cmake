set(AUTODIFF_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetAutodiff)
	find_path(AUTODIFF_DIR NAMES autodiff PATHS ${AUTODIFF_CMAKE_DIR}/../_deps/autodiff-src)

  if (NOT AUTODIFF_DIR ) 
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR "${AUTODIFF_CMAKE_DIR}/../_deps")
    FetchContent_Declare(autodiff
    GIT_REPOSITORY https://github.com/autodiff/autodiff.git
    GIT_TAG v0.6.12
    )

    FetchContent_GetProperties(autodiff)
    
    if(NOT autodiff_POPULATED)
      message(STATUS "Autodiff library not found locally, downloading it.")
      FetchContent_Populate(autodiff)
    endif()

    set(AUTODIFF_DIR ${autodiff_SOURCE_DIR})
    message(STATUS "Using downloaded Autodiff library at: ${AUTODIFF_DIR}")

  else ()

    message(STATUS "Autodiff Library found: ${AUTODIFF_DIR}")

  endif()

  include_directories(${AUTODIFF_DIR})

endfunction()