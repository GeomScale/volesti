function(GetBoost)

	find_path(BOOST_DIR NAMES Eigen PATHS ../../external)

	if (NOT BOOST_DIR) 

    set(BOOST_URL "https://dl.bintray.com/boostorg/release/1.71.0/source/boost_1_71_0.tar.bz2" CACHE STRING "Boost download URL")
    set(BOOST_URL_SHA256 "d73a8da01e8bf8c7eda40b4c84915071a8c8a0df4a6734537ddde4a8580524ee" CACHE STRING "Boost download URL SHA256 checksum")

    include(FetchContent)
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
