set(LP_SOLVE_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetLPSolve)
if(NOT APPLE)
  find_path(LP_SOLVE_DIR NAMES lpsolve.h PATHS ${LP_SOLVE_CMAKE_DIR}/../_deps/lpsolve-src)
endif(NOT APPLE)

  if (NOT LP_SOLVE_DIR)
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${LP_SOLVE_CMAKE_DIR}/../_deps")
      FetchContent_Declare(
          lpsolve
          URL https://downloads.sourceforge.net/project/lpsolve/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz 
          URL_HASH MD5=a829a8d9c60ff81dc72ff52363703886
      )

      FetchContent_GetProperties(lpsolve)

      if(NOT lpsolve_POPULATED)
          message(STATUS "lp_solve library not found locally, downloading it.")
          FetchContent_MakeAvailable(lpsolve)
      endif()

      set(LP_SOLVE_DIR "${lpsolve_SOURCE_DIR}")
      message(STATUS "Using downloaded lp_solve at: ${LP_SOLVE_DIR}")

  else()

      message(STATUS "lp_solve library found: ${LP_SOLVE_DIR}")

  endif()

  # Main lp_solve header (needed globally for lp_oracles)
  include_directories (BEFORE ${LP_SOLVE_DIR})

  add_library (lp_solve
  ${LP_SOLVE_DIR}/bfp/bfp_LUSOL/lp_LUSOL.c
  ${LP_SOLVE_DIR}/bfp/bfp_LUSOL/LUSOL/lusol.c
  ${LP_SOLVE_DIR}/colamd/colamd.c
  ${LP_SOLVE_DIR}/ini.c
  ${LP_SOLVE_DIR}/shared/commonlib.c
  ${LP_SOLVE_DIR}/shared/mmio.c
  ${LP_SOLVE_DIR}/shared/myblas.c
  ${LP_SOLVE_DIR}/lp_crash.c
  ${LP_SOLVE_DIR}/lp_Hash.c
  ${LP_SOLVE_DIR}/lp_lib.c
  ${LP_SOLVE_DIR}/lp_matrix.c
  ${LP_SOLVE_DIR}/lp_MDO.c
  ${LP_SOLVE_DIR}/lp_mipbb.c
  ${LP_SOLVE_DIR}/lp_MPS.c
  ${LP_SOLVE_DIR}/lp_params.c
  ${LP_SOLVE_DIR}/lp_presolve.c
  ${LP_SOLVE_DIR}/lp_price.c
  ${LP_SOLVE_DIR}/lp_pricePSE.c
  ${LP_SOLVE_DIR}/lp_report.c
  ${LP_SOLVE_DIR}/lp_scale.c
  ${LP_SOLVE_DIR}/lp_simplex.c
  ${LP_SOLVE_DIR}/lp_SOS.c
  ${LP_SOLVE_DIR}/lp_utils.c
  ${LP_SOLVE_DIR}/lp_wlp.c)

  # Scoped compile definitions & includes — prevents leaking to SuiteSparse etc.
  target_compile_definitions(lp_solve PRIVATE
    YY_NEVER_INTERACTIVE
    LoadInverseLib=0
    LoadLanguageLib=0
    LoadableBlasLib=0
    RoleIsExternalInvEngine
    INVERSE_ACTIVE=3
  )
  target_include_directories(lp_solve PRIVATE
    ${LP_SOLVE_DIR}/bfp
    ${LP_SOLVE_DIR}/bfp/bfp_LUSOL
    ${LP_SOLVE_DIR}/bfp/bfp_LUSOL/LUSOL
    ${LP_SOLVE_DIR}/colamd
    ${LP_SOLVE_DIR}/shared
  )

endfunction()
