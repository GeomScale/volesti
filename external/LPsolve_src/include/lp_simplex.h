// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#ifndef HEADER_lp_simplex
#define HEADER_lp_simplex

#include "lp_types.h"

#define ForceDualSimplexInBB               /* Force use/switch of dual simplex in B&B */
#define AssumeHighAccuracyInBB    /* No iteration of simplex solves at infeasibility */
/*#define UseLongStepPruning*/
/*#define UseLongStepDualPhase1*/
#define primal_UseRejectionList
#define dual_UseRejectionList
#define dual_RemoveBasicFixedVars
/*#define dual_Phase1PriceEqualities */   /* Force elimination of equality slacks */
#define AcceptMarginalAccuracy

#ifdef __cplusplus
extern "C" {
#endif

/* Put function headers here */
STATIC int primloop(lprec *lp, MYBOOL primalfeasible, LPSREAL primaloffset);
STATIC int dualloop(lprec *lp, MYBOOL dualfeasible, int dualinfeasibles[], LPSREAL dualoffset);
STATIC int spx_run(lprec *lp, MYBOOL validInvB);
STATIC int spx_solve(lprec *lp);
STATIC int lag_solve(lprec *lp, LPSREAL start_bound, int num_iter);
STATIC int heuristics(lprec *lp, int mode);
STATIC int lin_solve(lprec *lp);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_simplex */

