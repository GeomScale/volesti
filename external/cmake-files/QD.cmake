set(QD_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetQD)
  find_path(QD_DIR NAMES config.h PATHS ${QD_CMAKE_DIR}/../_deps/qd-src/)

  if (NOT QD_DIR)
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${QD_CMAKE_DIR}/../_deps")
      FetchContent_Declare(
          qd
          URL https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz
      )

      FetchContent_GetProperties(qd)

      if(NOT qd_POPULATED)
          message(STATUS "QD library not found locally, downloading it.")
          FetchContent_Populate(qd)
      endif()

      set(QD_DIR "${qd_SOURCE_DIR}")
      message(STATUS "Using downloaded QD at: ${QD_DIR}")

  else()

      message(STATUS "QD library found: ${QD_DIR}")

  endif()

  include_directories(BEFORE "${QD_DIR}/include/")
  message(STATUS "configuring the QD library")
  execute_process(
    COMMAND ./configure
    WORKING_DIRECTORY ${QD_DIR}
    OUTPUT_FILE CMD_OUTPUT
    RESULT_VARIABLE EXECUTE
  )
  if(NOT ${EXECUTE} EQUAL "0")
    message(FATAL_ERROR "./configure QD library failed")
  endif()

  execute_process(
      COMMAND make
      WORKING_DIRECTORY ${QD_DIR}
      OUTPUT_FILE qd_compilation.txt
      RESULT_VARIABLE EXECUTE_MAKE
  )

  if(NOT ${EXECUTE_MAKE} EQUAL "0")
    message(FATAL_ERROR "building the QD library failed")
  endif()

  find_library(QD_LIB NAMES libqd.a PATHS "${QD_DIR}/src/.libs")

endfunction()
