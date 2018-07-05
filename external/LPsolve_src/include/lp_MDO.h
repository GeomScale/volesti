// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#ifndef HEADER_MDO
#define HEADER_MDO

#include "lp_types.h"


#ifdef __cplusplus
extern "C" {
#endif

int __WINAPI getMDO(lprec *lp, MYBOOL *usedpos, int *colorder, int *size, MYBOOL symmetric);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_MDO */

