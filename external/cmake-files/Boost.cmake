set(BOOST_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetBoost)
	find_path(BOOST_DIR NAMES boost PATHS ${BOOST_CMAKE_DIR}/../_deps/boost-src)

	if (NOT BOOST_DIR) 
  
    set(BOOST_URL "https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2" CACHE STRING "Boost download URL")
    set(BOOST_URL_SHA256 "f0397ba6e982c4450f27bf32a2a83292aba035b827a5623a14636ea583318c41" CACHE STRING "Boost download URL SHA256 checksum")

    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR "${BOOST_CMAKE_DIR}/../_deps")
    FetchContent_Declare(
      Boost
      URL ${BOOST_URL}
      URL_HASH SHA256=${BOOST_URL_SHA256}
    )
    FetchContent_GetProperties(Boost)

    if(NOT Boost_POPULATED)
      message(STATUS "Fetching Boost")
      FetchContent_Populate(Boost)
      message(STATUS "Fetching Boost - done")
      set(BOOST_DIR ${boost_SOURCE_DIR})
    endif()

    message(STATUS "Using downloaded Boost library at ${BOOST_DIR}")

  else ()
    message(STATUS "Boost Library found: ${BOOST_DIR}")

  endif()

  include_directories(${BOOST_DIR})
  include_directories(${BOOST_DIR}/boost)

endfunction()
