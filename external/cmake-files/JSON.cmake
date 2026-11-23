set(JSON_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetJSON)
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR "${JSON_CMAKE_DIR}/../_deps")
    FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )
    message(STATUS "Fetching JSON")
    FetchContent_MakeAvailable(json)
endfunction()