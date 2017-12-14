/*
FILE DETAILS
description: header file needed for matrix operations in LCP solver, and some macros for the main
LCP function 
project: MPT 3.0
filename: lcp_matrix.h
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
details: Included BLAS and LAPACK libraries which are standard in MATLAB. Replaced int 
declarations with ptrdiff_t to be 32/64 compatible. Added declaration of LCP_options structure. 
*/


#ifndef LCP_MATRIX__H
#define LCP_MATRIX__H

/**** MODULES ****/
/* Standard headers */
#include <stdlib.h> 
#include <stdio.h>
#include <string.h> 
#include <stddef.h> /* to have ptrdiff_t */
#include <math.h>
#include <time.h>

/* BLAS and LAPACK headers are in $MATLABROOT/extern/include 
and are on the default search path when compiling with standard options */
#include <blas.h>
#include <lapack.h>


/**** MACROS ****/

/* Select an element from a column-major matrix */
#define C_SEL(mat,r,c) ((mat)->A[(r)+(c)*((mat)->rows_alloc)])

/* Select an element from a row-major matrix */
#define R_SEL(mat,r,c) ((mat)->A[(c)+(r)*((mat)->cols_alloc)])

/* Get a matrix */
#define MAT(mat) ((mat).A)

/* Get a matrix pointer */
#define pMAT(mat) ((mat)->A)

/* Get a pointer to a given column in a matrix */
#define COL(mat,n) &([C_SEL(mat,0,n)])

/* Get number of rows of a given matrix */
#define Matrix_Rows(pA) ((pA)->rows)

/* Get number of columns of a given matrix */
#define Matrix_Cols(pA) ((pA)->cols)

/* Dimension of an index set */
#define Index_Length(pI) ((pI)->len)

/* Assign value "x" to position "i" in the index set "I" */
#define Index_Set(pI,i,x) (((pI)->I[(i)]) = (x))

/* Obtain value of the index set "I" at the position "i" */
#define Index_Get(pI,i) ((pI)->I[(i)])


/**** TYPEDEFS AND STRUCTURES ****/

/* Matrices */
/* If a matrix is column-major, then the form is: */
/*      [M(0) M(rows_alloc)   M(2rows_alloc)   ... ] */
/*      [M(1) M(rows_alloc+1) M(2rows_alloc+1) ... ] */
/*      [...  ...             ....             ... ] */
typedef struct LCP_Matrix {
	double *A;
	ptrdiff_t rows,cols;
	ptrdiff_t rows_alloc, cols_alloc;
	char *name;
} T_Matrix, *PT_Matrix;


/* Index sets for basis defintion */
typedef struct LCP_IndexSet {
	/* C-based indexing used throughout. i.e. 0 is the first element */
	ptrdiff_t *I;
	ptrdiff_t len, len_alloc;
	char *name;
} T_Index, *PT_Index;

/* The basis */
typedef struct LCP_Basis {
	PT_Index pW;  /* Slacks in the basis */
	PT_Index pZ;  /* Variables in the basis */
	PT_Index pWc; /* Complement of W */

	PT_Index pP;  /* Work index (used to list positive variables) */
  
	/* B      = [I G;0 X]; */
	/* inv(B) = [I -G*iX; 0 iX] */
	/* G = M_W,Z */
	/* X = M_Wc,X */

	PT_Matrix pG, pX, pF, pLX, pUX, pUt;  /* F is used as factorization of X on exit from dgesv routine,
						                         Ut is full, column major upper triangular matrix */
	double *y,*z,*w, *v, *r;    /* work vectors */
	ptrdiff_t *p; /* work vector for dgesv, dgetrs */
	double *s; /* work vector for dgels */
} T_Basis, *PT_Basis;


/* options structure */
typedef struct LCP_options {
	double zerotol; /* zero tolerance*/
	double lextol; /* lexicographic tolerance - a small treshold to determine if values are equal */
	ptrdiff_t maxpiv; /* maximum number of pivots */
	ptrdiff_t nstepf; /* after n steps refactorize the basis */
	int clock;  /* print computation time */
	int verbose; /* verbose mode */
	int routine; /* routine from LAPACK library to solve linear equations in Basis_solve*/
	double timelimit; /* time limit for iterations in seconds */
	int normalize; /* normalize matrices M, q */
    double normalizethres; /* normalization threshold */
} T_Options;


/**** FUNCTIONS ****/

/* Basis functions */
PT_Basis Basis_Init(ptrdiff_t m);
void Basis_Free(PT_Basis pB);
void Basis_Print(PT_Basis pB);
ptrdiff_t  Basis_Pivot(PT_Basis pB, PT_Matrix pA, ptrdiff_t enter, ptrdiff_t leave);
ptrdiff_t Basis_Solve(PT_Basis pB, double *q, double *x, T_Options options);
void Reinitialize_Basis(ptrdiff_t m, PT_Basis pB);

void LU_Expand(PT_Matrix L, PT_Matrix U, double *y, double *z, double *w);
void LU_Replace_Col(PT_Matrix L, PT_Matrix U, ptrdiff_t kcol, double *y, double *z, double *w);
void LU_Replace_Row(PT_Matrix L, PT_Matrix U, ptrdiff_t krow, double *y, double *z, double *w);
void LU_Shrink(PT_Matrix L, PT_Matrix U, ptrdiff_t krow, ptrdiff_t kcol, double *y, double *z, double *w);
void LU_Solve0(PT_Matrix pL, PT_Matrix pU, double *vec, double *x);
void LU_Solve1(PT_Matrix pL, PT_Matrix pU, double *vec, double *x, ptrdiff_t *info);
void LU_Refactorize(PT_Basis pB);

/* Matrix functions */
PT_Matrix Matrix_Init(ptrdiff_t m, ptrdiff_t n, const char *name);
PT_Matrix Matrix_Init_NoAlloc(ptrdiff_t m, ptrdiff_t n, const char *name);
void Matrix_Free(PT_Matrix pMat);
void Matrix_Free_NoAlloc(PT_Matrix pMat);
void Matrix_Print_col(PT_Matrix pMat);
void Matrix_Print_row(PT_Matrix pMat);
void Matrix_Print_utril_row(PT_Matrix pMat);
void Matrix_Print_size(PT_Matrix pMat);
void Matrix_AddRow(PT_Matrix pA, double *row);
void Matrix_AddCol(PT_Matrix pA, double *col);
void Matrix_RemoveRow(PT_Matrix pA, ptrdiff_t row);
void Matrix_RemoveCol(PT_Matrix pA, ptrdiff_t col);
void Matrix_GetRow(PT_Matrix pA, double *row, ptrdiff_t row_index, PT_Index col_index);
void Matrix_GetCol(PT_Matrix pA, double *col, PT_Index row_index, ptrdiff_t col_index);
void Matrix_ReplaceCol(PT_Matrix pA, ptrdiff_t kcol, double *vec);
void Matrix_ReplaceRow(PT_Matrix pA, ptrdiff_t krow, double *vec);
void Matrix_Copy(PT_Matrix pA, PT_Matrix pB, double *w);
void Matrix_Triu(PT_Matrix pF, PT_Matrix pU, double *w);
void Matrix_Uzeros(PT_Matrix pF);
void Matrix_Full2RowTriangular(PT_Matrix pF, PT_Matrix  pU, double *w);
void Matrix_Transpose(PT_Matrix pA, PT_Matrix pB, double *w);
void Vector_Print_raw(double *vec, ptrdiff_t m);

/* Index set functions */
PT_Index Index_Init(ptrdiff_t len, const char *name);
void Index_Free(PT_Index pIndex);
void Index_Print(PT_Index pI);
void Index_Remove(PT_Index pI, ptrdiff_t leave);
void Index_Add(PT_Index pI, ptrdiff_t num);
void Index_Replace(PT_Index pI, ptrdiff_t index, ptrdiff_t newval);
ptrdiff_t Index_FindElement(PT_Index pI, ptrdiff_t element);

#endif
