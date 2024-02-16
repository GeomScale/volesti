set(BOOST_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetBoost)
	find_path(BOOST_DIR NAMES boost PATHS ${BOOST_CMAKE_DIR}/../_deps/boost-src)

	if (NOT BOOST_DIR)

    set(BOOST_URL "https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.bz2" CACHE STRING "Boost download URL")
    set(BOOST_URL_SHA256 "cc4b893acf645c9d4b698e9a0f08ca8846aa5d6c68275c14c3e7949c24109454" CACHE STRING "Boost download URL SHA256 checksum")

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
