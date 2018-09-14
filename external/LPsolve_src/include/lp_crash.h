// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#ifndef HEADER_lp_crash
#define HEADER_lp_crash


#include "lp_types.h"

#define CRASH_SIMPLESCALE       /* Specify if we should use a simple absolute scaling threshold */

#define CRASH_THRESHOLD  0.167
#define CRASH_SPACER        10
#define CRASH_WEIGHT     0.500



#ifdef __cplusplus
__EXTERN_C {
#endif

STATIC MYBOOL crash_basis(lprec *lp);

#ifdef __cplusplus
}
#endif

#endif /* HEADER_lp_crash */

