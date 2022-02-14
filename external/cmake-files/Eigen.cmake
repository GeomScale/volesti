set(EIGEN_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetEigen)
  find_path(EIGEN_DIR NAMES Eigen PATHS ${EIGEN_CMAKE_DIR}/../_deps/eigen-src)

  if (NOT EIGEN_DIR) 
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR "${EIGEN_CMAKE_DIR}/../_deps")
    FetchContent_Declare(
      eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG 3.4.0
    )

    FetchContent_GetProperties(eigen)
    
    if(NOT eigen_POPULATED)
      message(STATUS "Eigen library not found locally, downloading it.")
      FetchContent_Populate(eigen)
    endif()

    set(EIGEN_DIR ${eigen_SOURCE_DIR})
    message(STATUS "Using downloaded Eigen library at: ${EIGEN_DIR}")

  else ()

    message(STATUS "Eigen Library found: ${EIGEN_DIR}")

  endif()

  include_directories(${EIGEN_DIR})

endfunction()
