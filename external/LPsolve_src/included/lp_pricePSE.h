// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#ifndef HEADER_lp_pricePSE
#define HEADER_lp_pricePSE

#include "lp_types.h"

#define ApplySteepestEdgeMinimum

#ifdef __cplusplus
extern "C" {
#endif

/* Price norm management routines */
STATIC MYBOOL initPricer(lprec *lp);
INLINE MYBOOL applyPricer(lprec *lp);
STATIC void simplexPricer(lprec *lp, MYBOOL isdual);
STATIC void freePricer(lprec *lp);
STATIC MYBOOL resizePricer(lprec *lp);
STATIC LPSREAL getPricer(lprec *lp, int item, MYBOOL isdual);
STATIC MYBOOL restartPricer(lprec *lp, MYBOOL isdual);
STATIC MYBOOL updatePricer(lprec *lp, int rownr, int colnr, LPSREAL *pcol, LPSREAL *prow, int *nzprow);
STATIC MYBOOL verifyPricer(lprec *lp);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_pricePSE */

