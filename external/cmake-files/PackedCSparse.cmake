set(PackedCSparse_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetPackedCSparse)
  find_path(PackedCSparse_DIR NAMES PackedCSparse PATHS "${PackedCSparse_CMAKE_DIR}/../_deps/polytopesampler-src/code/solver")
  if (NOT PackedCSparse_DIR)
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${PackedCSparse_CMAKE_DIR}/../_deps")
      FetchContent_Declare(
          polytopesampler
          GIT_REPOSITORY https://github.com/ConstrainedSampler/PolytopeSamplerMatlab.git
          GIT_TAG f3aaacef7a2aa7e906d7b73e1605845bbdf6d35b
      )
      FetchContent_GetProperties(polytopesampler)

      if(NOT polytopesampler_POPULATED)
          message(STATUS "PackedCSparse library not found locally, downloading it.")
          FetchContent_Populate(polytopesampler)
      endif()

      set(PackedCSparse_DIR "${polytopesampler_SOURCE_DIR}/code/solver")
      message(STATUS "Using downloaded PackedCSparse at: ${PackedCSparse_DIR}/PackedCSparse")

  else()
      message(STATUS "PackedCSparse library found: ${PackedCSparse_DIR}")

  endif()

  include_directories(BEFORE "${PackedCSparse_DIR}")

  set(QDSOURCES
  ${PackedCSparse_DIR}/qd/bits.cc
  ${PackedCSparse_DIR}/qd/c_dd.cc
  ${PackedCSparse_DIR}/qd/c_qd.cc
  ${PackedCSparse_DIR}/qd/dd_const.cc
  ${PackedCSparse_DIR}/qd/dd_real.cc
  ${PackedCSparse_DIR}/qd/fpu.cc
  ${PackedCSparse_DIR}/qd/qd_const.cc
  ${PackedCSparse_DIR}/qd/qd_real.cc
  ${PackedCSparse_DIR}/qd/util.cc
  )
add_library(QD_LIB ${QDSOURCES})

endfunction()
