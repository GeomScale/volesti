set(QD_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetQD)
  set(QD_DIR "${QD_CMAKE_DIR}/../PackedCSparse/qd")
  include_directories(BEFORE "${QD_DIR}")
  set(QD_SOURCES
    ${QD_DIR}/bits.cc
    ${QD_DIR}/c_dd.cc
    ${QD_DIR}/c_qd.cc
    ${QD_DIR}/dd_const.cc
    ${QD_DIR}/dd_real.cc
    ${QD_DIR}/fpu.cc
    ${QD_DIR}/qd_const.cc
    ${QD_DIR}/qd_real.cc
    ${QD_DIR}/util.cc
        )
  add_library(QD_LIB ${QD_SOURCES})

endfunction()
