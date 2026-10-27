set(HIGHS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetHighs)
    # Prefer an installed HiGHS version.
    find_package(highs CONFIG QUIET)
    if(highs_FOUND)
        message(STATUS "Using installed HiGHS ${highs_VERSION}")
        return()
    endif()

    include(FetchContent)

    # Build highs without tests, examples or shared library.
    set(BUILD_TESTING OFF)
    set(BUILD_EXAMPLES OFF)
    set(BUILD_SHARED_LIBS OFF)

    FetchContent_Declare(
        highs
        GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
        GIT_TAG v1.10.0
        GIT_SHALLOW TRUE
    )
    message(STATUS "Fetching HiGHS")
    FetchContent_MakeAvailable(highs)
endfunction()