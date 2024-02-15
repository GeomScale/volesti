set(LP_SOLVE_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
function(GetLPSolve)
  find_path(LP_SOLVE_DIR NAMES lpsolve.h PATHS ${LP_SOLVE_CMAKE_DIR}/../_deps/lpsolve-src)

  if (NOT LP_SOLVE_DIR)
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${LP_SOLVE_CMAKE_DIR}/../_deps")
      FetchContent_Declare(
          lpsolve
          URL https://webwerks.dl.sourceforge.net/project/lpsolve/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz
          URL_HASH MD5=a829a8d9c60ff81dc72ff52363703886
      )

      FetchContent_GetProperties(lpsolve)

      if(NOT lpsolve_POPULATED)
          message(STATUS "lp_solve library not found locally, downloading it.")
          FetchContent_Populate(lpsolve)
      endif()

      set(LP_SOLVE_DIR "${lpsolve_SOURCE_DIR}")
      message(STATUS "Using downloaded lp_solve at: ${LP_SOLVE_DIR}")

  else()

      message(STATUS "lp_solve library found: ${LP_SOLVE_DIR}")

  endif()

  #to disable interactive mode of lp_solve lex parser
  add_compile_options(-DYY_NEVER_INTERACTIVE)

  add_compile_options(-DLoadInverseLib=0)
  add_compile_options(-DLoadLanguageLib=0)
  add_compile_definitions(RoleIsExternalInvEngine)
  add_compile_definitions(INVERSE_ACTIVE=3)

  include_directories (BEFORE ${LP_SOLVE_DIR})
  include_directories (BEFORE ${LP_SOLVE_DIR}/bfp)
  include_directories (BEFORE ${LP_SOLVE_DIR}/bfp/bfp_LUSOL)
  include_directories (BEFORE ${LP_SOLVE_DIR}/bfp/bfp_LUSOL/LUSOL)
  include_directories (BEFORE ${LP_SOLVE_DIR}/colamd)
  include_directories (BEFORE ${LP_SOLVE_DIR}/shared)

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

endfunction()
