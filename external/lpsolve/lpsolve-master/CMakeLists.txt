# Copyright (C) 2014 Tai Chi Minh Ralph Eastwood.
# Redistribution and use of this file is allowed according to the terms of the MIT license.
# For details see the COPYRIGHT file distributed with Flymir LuaDist.
# Please note that the package source code is licensed under its own license.

project(lpsolve C)
cmake_minimum_required(VERSION 2.8)

include (cmake/dist.cmake)
include (cmake/lpsolve.cmake)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/shared)
include_directories(${PROJECT_SOURCE_DIR}/bfp)
include_directories(${PROJECT_SOURCE_DIR}/bfp/bfp_LUSOL)
include_directories(${PROJECT_SOURCE_DIR}/bfp/bfp_LUSOL/LUSOL)
include_directories(${PROJECT_SOURCE_DIR}/colamd)

add_definitions(-DYY_NEVER_INTERACTIVE -DPARSER_LP -DINVERSE_ACTIVE=INVERSE_LUSOL -DRoleIsExternalInvEngine)
add_definitions(-DLLONG=long)

set(LPSOLVE_SOURCES
	lp_MDO.c
	shared/commonlib.c
	shared/mmio.c
	shared/myblas.c
	ini.c
	fortify.c
	colamd/colamd.c
	lp_rlp.c
	lp_crash.c
	bfp/bfp_LUSOL/lp_LUSOL.c
	bfp/bfp_LUSOL/LUSOL/lusol.c
	lp_Hash.c
	lp_lib.c
	lp_wlp.c
	lp_matrix.c
	lp_mipbb.c
	lp_MPS.c
	lp_params.c
	lp_presolve.c
	lp_price.c
	lp_pricePSE.c
	lp_report.c
	lp_scale.c
	lp_simplex.c
	lp_SOS.c
	lp_utils.c
	yacc_read.c
)

set(LPSOLVE_HEADERS
	lp_bit.h
	lp_crash.h
	lp_explicit.h
	lp_fortify.h
	lp_Hash.h
	lp_lib.h
	lp_matrix.h
	lp_MDO.h
	lp_mipbb.h
	lp_MPS.h
	lp_presolve.h
	lp_price.h
	lp_pricePSE.h
	lp_report.h
	lp_rlp.h
	lp_scale.h
	lp_simplex.h
	lp_solveDLL.h
	lp_SOS.h
	lp_types.h
	lp_utils.h
	lp_wlp.h
)

add_library(lpsolve55 ${LPSOLVE_SOURCES})
target_link_libraries(lpsolve55 m dl)

add_executable(lp_solve ${LPSOLVE_SOURCES} lp_solve/lp_solve.c)
target_link_libraries(lp_solve m dl)

install(TARGETS lpsolve55 DESTINATION lib)
install(FILES ${LPSOLVE_HEADERS} DESTINATION include)
