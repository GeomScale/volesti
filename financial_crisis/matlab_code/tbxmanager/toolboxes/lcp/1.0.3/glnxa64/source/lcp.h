/*
FILE DETAILS
description: headers and declarations to the main LCP function
project: MPT 3.0
filename: lcp.h
author: Colin N. Jones, Automatic Control Laboratory, ETH Zurich, 2006
 
LICENSE:
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

REVISION HISTORY:
date: Nov, 2010
revised by: Martin Herceg, Automatic Control Laboratory, ETH Zurich, 2010
details: Removed definition of ZEROTOL which will be passed in options structure by the user.
*/

#ifndef LCP_MAIN__H
#define LCP_MAIN__H

/*** MODULES ***/
/* Matlab mex interface */
#ifdef  MATLAB_MEX_FILE
#include "mex.h"
#endif

/* LCP functions */
#include "lcp_matrix.h"

/**** DEFINITIONS AND MACROS ****/

#define LCP_FEASIBLE    1
#define LCP_INFEASIBLE -1
#define LCP_UNBOUNDED  -2
#define LCP_PRETERMINATED  -3
#define LCP_CODEERROR  -4

/* #define ASSERT(t) t?:mexErrMsgTxt("Assertion failed"); */


/*** FUNCTIONS ***/

PT_Matrix lcp_Matrix_Init(ptrdiff_t m, double *M, double *q);
int lcp(PT_Basis pB, PT_Matrix pA, ptrdiff_t *pivs, T_Options options);
void NormalizeMatrix (ptrdiff_t m, ptrdiff_t n, double *A, double *b, double *r, double *c,  T_Options options);

#endif
