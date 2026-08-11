set(HIGHS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetHighs)
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR "${HIGHS_CMAKE_DIR}/../_deps")
    FetchContent_Declare(
        highs
        GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
        GIT_TAG v1.10.0
    )
    message(STATUS "Fetching HiGHS")
    FetchContent_MakeAvailable(highs)
endfunction()