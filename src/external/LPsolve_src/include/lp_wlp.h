// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#ifndef HEADER_lp_lp
#define HEADER_lp_lp

#include "lp_types.h"


#ifdef __cplusplus
extern "C" {
#endif

/* Put function headers here */
MYBOOL LP_writefile(lprec *lp, char *filename);
MYBOOL LP_writehandle(lprec *lp, FILE *output);


#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_lp */

