set(SPECTRA_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})

function(GetSpectra)
  # Use the bundled Spectra library included in external/Spectra/
  set(BUNDLED_SPECTRA_DIR "${SPECTRA_CMAKE_DIR}/../Spectra/include")
  if (EXISTS "${BUNDLED_SPECTRA_DIR}/Spectra/SymEigsSolver.h")
    set(SPECTRA_DIR ${BUNDLED_SPECTRA_DIR})
    message(STATUS "Using bundled Spectra library at: ${SPECTRA_DIR}")
  else()
    # Fall back to downloading if bundled version is missing
    find_path(SPECTRA_DIR NAMES Spectra/SymEigsSolver.h PATHS ${SPECTRA_CMAKE_DIR}/../_deps/spectra-src PATH_SUFFIXES include)
    if (NOT SPECTRA_DIR)
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${SPECTRA_CMAKE_DIR}/../_deps")
      FetchContent_Declare(
        spectra
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG master
      )
      FetchContent_GetProperties(spectra)
      if(NOT spectra_POPULATED)
        message(STATUS "Spectra library not found locally, downloading it.")
        FetchContent_MakeAvailable(spectra)
      endif()
      set(SPECTRA_DIR ${spectra_SOURCE_DIR}/include)
      message(STATUS "Using downloaded Spectra library at: ${SPECTRA_DIR}")
    else ()
      message(STATUS "Spectra Library found: ${SPECTRA_DIR}")
    endif()
  endif()
  include_directories(${SPECTRA_DIR})
endfunction()