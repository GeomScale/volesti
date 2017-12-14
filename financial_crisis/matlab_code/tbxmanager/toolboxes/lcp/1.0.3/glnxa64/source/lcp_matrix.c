/*
FILE DETAILS
description: help routines for matrix operations, index set operations,
and some macros for the main LCP function
project: MPT 3.0
filename: lcp_matrix.c
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
date: Aug, 2011
revised by: Martin Herceg
details: Due to problems in freeing vectors, we allocate vectors in basis 
          always 5 elements more that it is needed.
 
date: Nov, 2010
revised by: Martin Herceg, Automatic Control Laboratory, ETH Zurich, 2010
details: Basic solution involved ih the Basis_solve function is performed via LUMOD routine (by
default) and 2 routines from LAPACK library. User can choose the routine depending on the
options. Added testing on singularity. If the LUMOD factorization is not satisfactory (i.e. matrix
might be singular), L,U factors are recomputed with dgetrf routine from LAPACK.

*/


/************************************************************************
 MODULES 
************************************************************************/

/* standard C headers are included in "lcp_matrix.h" */
#include "lcp_matrix.h"
#include "lumod_dense.h"

/************************************************************************
 Basis functions
************************************************************************/

/* Initialize basis */
PT_Basis Basis_Init(ptrdiff_t m)
{
	ptrdiff_t i;
	PT_Basis pB;

	pB = (PT_Basis)malloc(sizeof(T_Basis));
  
	/* Initialize indices */
	pB->pW  = Index_Init(m,"W");
	pB->pZ  = Index_Init(m,"Z");
	pB->pWc = Index_Init(m,"Wc");
	pB->pP  = Index_Init(m,"P");

	/* Initial basis is I = [1:m] */
	/*   or W = [1:m], Z = [], Wc = [] */
	for(i=0;i<m;i++) (pB->pW)->I[i] = i;
	(pB->pZ)->len = 0;
	(pB->pWc)->len = 0;

	/* Reserve space for a G matrix of max size (m x m) */
	pB->pG = Matrix_Init(m, m, "B.G");
	/* G is of size |W|x|Z| */
	(pB->pG)->rows = m;
	(pB->pG)->cols = 0;

	/* Reserve space for an X matrix of max size (m x m) */ 
	/* this matrix is to be factored */
	/* X = inv(LX)*UX */
	pB->pX = Matrix_Init(m, m, "B.X");
	(pB->pX)->rows = pB->pWc->len;
	(pB->pX)->cols = pB->pZ->len;

	/* Reserve space for factored matrix F = P*L*U of max size (m x m) */ 
	pB->pF = Matrix_Init(m, m, "B.F");
	(pB->pF)->rows = pB->pWc->len;
	(pB->pF)->cols = pB->pZ->len;

	/* Setup the initial decomposition of X for LUmod */
	/* X starts empty - so just reserve space */
	/* LUmod performs factorization of type L*X = U */
	/* LX is row major, full */
	/* UX is row major, upper triangular */
	/* Ut is column major, full with upper triangular entries */
	pB->pLX = Matrix_Init(m,m,"B.LX");
	pB->pUX = Matrix_Init(m,m,"B.UX");
	pB->pUt = Matrix_Init(m,m,"B.Ut");
	pB->pLX->rows = pB->pWc->len;
	pB->pLX->cols = pB->pZ->len;
	pB->pUX->rows = pB->pWc->len;
	pB->pUX->cols = pB->pZ->len;
	pB->pUt->rows = pB->pWc->len;
	pB->pUt->cols = pB->pZ->len;
	

	/* Work vectors*/
    /* allocate always 5 elements more because there were problems detected
     * when freeing vectors */
	pB->y = (double *)calloc(m+5,sizeof(double));
	pB->z = (double *)calloc(m+5,sizeof(double));
	pB->w = (double *)calloc(m+5,sizeof(double));
	pB->v = (double *)calloc(m+5,sizeof(double));
	pB->r = (double *)calloc(m+5,sizeof(double));
	pB->p = (ptrdiff_t *)calloc(m+5,sizeof(ptrdiff_t)); /* for dgesv, dgetrs routine */
	pB->s = (double *)calloc(2*m+5,sizeof(double)); /* for dgels routine */

	return pB;
}

/* Free basis */
void Basis_Free(PT_Basis pB)
{
	Index_Free(pB->pW);
	Index_Free(pB->pZ);
	Index_Free(pB->pWc);
	Index_Free(pB->pP);

	Matrix_Free(pB->pG);
	Matrix_Free(pB->pX);
	Matrix_Free(pB->pF);
	Matrix_Free(pB->pLX);
	Matrix_Free(pB->pUX);
	Matrix_Free(pB->pUt);
	
	if (pB->pW!=NULL)
		free(pB->pW);
	if (pB->pZ!=NULL)
		free(pB->pZ);
	if (pB->pWc!=NULL)
		free(pB->pWc);
	if (pB->pP!=NULL)
		free(pB->pP);
	if (pB->pG!=NULL)
		free(pB->pG);
	if (pB->pX!=NULL)
		free(pB->pX);
	if (pB->pF!=NULL)
		free(pB->pF);
	if (pB->pLX!=NULL)
		free(pB->pLX);
	if (pB->pUX!=NULL)
		free(pB->pUX);
	if (pB->pUt!=NULL)
		free(pB->pUt);

	if (pB->y!=NULL)
		free(pB->y);
	if (pB->z!=NULL)
		free(pB->z);
	if (pB->w!=NULL)
		free(pB->w);
	if (pB->v!=NULL)
		free(pB->v);
	if (pB->r!=NULL)
		free(pB->r);
	if (pB->p!=NULL)
		free(pB->p);
	if (pB->s!=NULL)
		free(pB->s);
	
}

/* Print basis */
void Basis_Print(PT_Basis pB)
{
	printf("Current basis:\n");
	Index_Print(pB->pW);
	Index_Print(pB->pZ);
	Index_Print(pB->pWc);
  
	Matrix_Print_col(pB->pG);
	Matrix_Print_col(pB->pX);
	/* Matrix_Print_col(pB->pF); */
	Matrix_Print_utril_row(pB->pUX);
	Matrix_Print_row(pB->pLX);

}

/* The main function: pivot */
ptrdiff_t  Basis_Pivot(PT_Basis pB, PT_Matrix pA, ptrdiff_t enter, ptrdiff_t leave)
{
	ptrdiff_t m, nW, Wl, c, r;
	PT_Index pW, pZ, pWc;
	PT_Matrix pG, pX, pF;
	ptrdiff_t left;

	/* abbreviations */
	pW  = pB->pW;
	pWc = pB->pWc;
	pZ  = pB->pZ;
	pG  = pB->pG;
	pX = pB->pX;
	pF = pB->pF;

	/* get dimensions */
	m   = Matrix_Rows(pA);
	nW  = Index_Length(pW);

	/* determine the leaving variable */
	if(leave < Index_Length(pW))
	{
		left = Index_Get(pW,leave);
	}
	else
	{
		left = Index_Get(pZ,leave-Index_Length(pW))+m;
	}
  
/*   printf("\n ---------------------------------------\n"); */
/*   printf("enter = %ld, leave = %ld\n", enter, leave); */
  
	/* Four cases */
	if(enter < m)
	{
		/* r = find(B.Wc == enter) */
		r = Index_FindElement(pWc, enter);

		if(leave < nW)
		{
			/* printf("case 1:\n"); */
			Wl = Index_Get(pW,leave);

			/* G_l,* = M_e,Z */
			/* Basis_Print(pB); */
			/* Update matrix G with the row corresponding to the entering variable */
			Matrix_GetRow(pA, pB->y, enter, pZ);
			Matrix_ReplaceRow(pG, leave, pB->y);

			/* X_r,* = M_Wl,Z */
			/* Update matrix X with the row corresponing to the leaving variable 
			   from index set pW. Similarly, we must update matrix F as its gets overwritten
			   by dgesv routine.   */
			Matrix_GetRow(pA, pB->y, Wl, pZ);
			Matrix_ReplaceRow(pF, r, pB->y);      
			Matrix_ReplaceRow(pX, r, pB->y);
			/* smart update with LUmod */
			LU_Replace_Row(pB->pLX, pB->pUX, r, pB->y, pB->z, pB->w);

			/* Wc_r = Wl, W_l = e */
			/* update index sets */
			Index_Replace(pWc, r, Wl);
			Index_Replace(pW,  leave, enter);

		}
		else /* leave >= nW */
		{
			/* printf("case 2:\n"); */
			Index_Add(pW,enter);
			c = leave-nW;

			/* Replace row r and column c with the last row and column of X
			   and remove last row/column. */
			Matrix_RemoveRow(pX,r);
			Matrix_RemoveCol(pX,c);
			/* We must update matrix F as well */
			Matrix_RemoveRow(pF,r);
			Matrix_RemoveCol(pF,c);
			/* smart update with LUmod */
			LU_Shrink(pB->pLX, pB->pUX, r, c, pB->y, pB->z, pB->w);

			/* Update indices */
			Index_Remove(pWc, r);
			Index_Remove(pZ,  c);

			/* Update G matrix */
			Matrix_RemoveCol(pG, c);
			Matrix_GetRow(pA, pB->y, enter, pZ);
			Matrix_AddRow(pG, pB->y);
		}
	}
	else /* enter >= m */
	{
		if(leave < nW)
		{
			/* printf("case 3:\n"); */
			Wl = Index_Get(pW,leave);

			/* Add row Wl and column e-m to X */
			if(Index_Length(pZ) > 0)
			{
				Matrix_GetRow(pA, pB->y, Wl, pZ);
				Matrix_GetCol(pA, pB->z, pWc, enter-m);
			}
			pB->y[Index_Length(pZ)] = C_SEL(pA,Wl,enter-m);

			/* updating F */
			Matrix_AddCol(pF, pB->z);
			Matrix_AddRow(pF, pB->y);

			/* updating X */
			Matrix_AddCol(pX, pB->z);
			Matrix_AddRow(pX, pB->y);
			/* update with LUmod */
			LU_Expand(pB->pLX, pB->pUX, pB->y, pB->z, pB->w);

			/* printf("pX col=\n"); */
			/* Matrix_Print_col(pX); */
     
			/* Remove row Wl from G */
			/* and add column e-m */
			Matrix_RemoveRow(pG, leave);
			Index_Remove(pW, leave);

			Matrix_GetCol(pA, pB->z, pW, enter-m);
			Matrix_AddCol(pG, pB->z);

			/* Extend Wc and Z */
			Index_Add(pWc, Wl);
			Index_Add(pZ, enter-m);

		}
		else /* leave >= nW */
		{
			/* printf("case 4:\n"); */
			/* Replace column */
			c = leave - nW;

			/* G_*,c = M_w,(e-m) */
			Matrix_GetCol(pA, pB->y, pW, enter-m);
			Matrix_ReplaceCol(pG, c, pB->y);

			/* X_*,c = M_Wc,(e-m) */
			Matrix_GetCol(pA, pB->z, pWc, enter-m);

			/* Vector_Print_raw(pB->z,m); */
			/* Basis_Print(pB); */

			/* update matrices X, F by replacing column */
			Matrix_ReplaceCol(pX, c, pB->z);
			Matrix_ReplaceCol(pF, c, pB->z);
			/* update for LUmod */
			LU_Replace_Col(pB->pLX, pB->pUX, c, pB->y, pB->z, pB->w);

			/* Update Z */
			Index_Replace(pZ, c, enter-m);
		}
	}

	/* create column-wise upper triangular matrix Ut from U*/
	Matrix_Triu(pB->pUt, pB->pUX, pB->r);

	return left;
}

/* Compute x = inv(B)*q for some basis and q */
/* Work variable B.z is modified */
/* Return information if the solution was computed */
ptrdiff_t Basis_Solve(PT_Basis pB, double *q, double *x, T_Options options)
{
/*   Basis B */
/*     I = W union (Z+m) -> ordered */
/*     X = M_/W,Z        <- invertible */
/*     G = M_W,Z */
/*   Bx_B = q -> x = [w_W;z_Z] */
/*     z_Z = iX*q_/W */
/*     w_W = q_W - G*x_/W */

	ptrdiff_t i, incx=1, nb, m, diml, info;
	char T;
	double alpha=-1.0, beta=1.0;


	/* x(1:length(W)) = q(W) */
	for(i=0;i<Index_Length(pB->pW);i++) x[i] = q[Index_Get(pB->pW, i)];

	/* X = [] */
	if(Index_Length(pB->pWc) == 0)
		return 0;

	/* z = q(B.Wc) */
	for(i=0;i<Index_Length(pB->pWc);i++) pB->z[i] = q[Index_Get(pB->pWc, i)];

	/* because the right hand side will be overwritten in dgesv, store it into 
	   temporary vector r */
	dcopy(&(Index_Length(pB->pWc)), pB->z, &incx, pB->r, &incx);

	/* printf("x :\n"); */
	/* Vector_Print_raw(pB->z,Index_Length(pB->pW)); */

	/* printf("z :\n"); */
	/* Vector_Print_raw(pB->z,Index_Length(pB->pWc)); */
 
	/* x = [x1;x2] */
	/* x2 = inv(X)*q(Wc) = inv(X)*z */
	nb = 1; /* number of right hand side in A*x=b is 1 */
	m = Matrix_Rows(pB->pF);

	/* Find solution to general problem A*x=b using LAPACK
	   All the arguments have to  be pointers.  A and b will be altered on exit. */
	switch (options.routine) {
	case 1:
		/* DGESV method implements  LU factorization of A. */		
		dgesv(&m, &nb, pMAT(pB->pF), &((pB->pF)->rows_alloc), 
		      pB->p, pB->r, &m, &info);
		break;
	case 2: 
		/*  DGELS solves overdetermined or underdetermined real linear systems
		    involving an M-by-N matrix A, or its transpose, using a QR or LQ
		    factorization of A.  It is assumed that A has full rank. */
		T = 'n'; /* Not transposed */
		diml = 2*m;
		dgels(&T, &m, &m, &nb, pMAT(pB->pF), &((pB->pF)->rows_alloc),
		      pB->r, &m, pB->s, &diml, &info );
		break;
	default : 
		/* solve with factors F or U */

		/* solve using LUMOD (no checks) */
		/* LU_Solve0(pB->pLX, pB->pUX, pB->z, &(x[Index_Length(pB->pW)])); */

		/* solve using LAPACK (contains also singularity checks) */
		/* need to pass info as a pointer, otherwise weird numbers are assigned */
		LU_Solve1(pB->pLX, pB->pUt, pB->z, &(x[Index_Length(pB->pW)]), &info);

		/* if something went wrong, refactorize and solve again */
		if (info!=0) {
			/* printf("info=%ld,  refactoring\n", info); */

			/* Matrix_Print_row(pB->pLX); */
			/* Matrix_Print_utril_row(pB->pUX); */
			/* Matrix_Print_col(pB->pUt); */

			LU_Refactorize(pB);
			/* if this fails again, then no minimum ratio is found -> exit with
			 * numerical problems flag */
			LU_Solve1(pB->pLX, pB->pUt, pB->z, &(x[Index_Length(pB->pW)]), &info);			
		}

	}
	

	/* x1 = -G*x2 + q(W) */
	/*  y = alpha*G*x + beta*y */
	/* alpha = -1.0; */
	/* beta  =  1.0; */
	T     = 'n'; /* Not transposed */
	if (options.routine>0) {
		/* take LAPACK solution */

		/* printf("lapack solution z:\n"); */
		/* Vector_Print_raw(pB->r,Index_Length(pB->pWc)); */

		/* matrix F after solution is factored in [L\U], we want the original format for the next call
		   to dgesv, thus create a copy F <- X */
		Matrix_Copy(pB->pX, pB->pF, pB->w);
 
		/* printf("x after lapack solution:\n"); */
		/* Vector_Print_raw(x,Index_Length(pB->pW)+Index_Length(pB->pWc)); */

		/* recompute the remaining variables according to basic solution */
		dgemv(&T, &(Matrix_Rows(pB->pG)), &(Matrix_Cols(pB->pG)),
		      &alpha, pMAT(pB->pG), &((pB->pG)->rows_alloc),
		      pB->r, &incx, &beta, x, &incx);
		/* append the basic solution at the end of x vector */
		dcopy(&(Index_Length(pB->pWc)), pB->r, &incx, &(x[Index_Length(pB->pW)]), &incx);
	} else {
		/* take LUmod solution */
		dgemv(&T, &(Matrix_Rows(pB->pG)), &(Matrix_Cols(pB->pG)),
		      &alpha, pMAT(pB->pG), &((pB->pG)->rows_alloc),
		      &(x[Index_Length(pB->pW)]), &incx,
		      &beta, x, &incx);
	}
	/* printf("y:\n"); */
	/* Vector_Print_raw(x,Matrix_Rows(pB->pG)); */

	return info;

}


/* reinitialize basis for repetitive calling */
void Reinitialize_Basis(ptrdiff_t m, PT_Basis pB)
{
	ptrdiff_t i;
	
	/* reset dimensions of index sets */
	for (i=0; i<m; i++) {
		(pB->pW)->I[i] = i;
	}
	pB->pW->len = m;
	pB->pP->len = m;
	pB->pWc->len = 0;
	pB->pZ->len = 0;

	/* reset dimension of all matrices */
	pB->pG->rows = m;
	pB->pG->cols = 0;
     
	pB->pX->rows = 0;
	pB->pX->cols = 0;

	pB->pF->rows = 0;
	pB->pF->cols = 0;

	pB->pLX->rows = 0;
	pB->pLX->cols = 0;
     
	pB->pUX->rows = 0;
	pB->pUX->cols = 0;
     
	pB->pUt->rows = 0;
	pB->pUt->cols = 0;

}

/************************************************************************
 LU Operations 
************************************************************************/

/* Add row and col to X */
/* They're added to the last row and column */
/* y,z,w = work of size double(m) */
/*         On entry: */
/*               y(*), z(*)  contain the new row and column. */
/*               y(n)        contains the new diagonal element of C. */
/*         On exit: */
/*               y(*)        is altered. */
void LU_Expand(PT_Matrix L, PT_Matrix U, double *y, double *z, double *w)
{
	int mode;
	ptrdiff_t tmp=0; /* not used */

	/* Increase size of matrix */
	L->rows = L->rows + 1;
	L->cols = L->cols + 1;
	U->rows = U->rows + 1;
	U->cols = U->cols + 1;

/*   ASSERT(L->rows <= L->rows_alloc); */
/*   ASSERT(U->rows <= U->rows_alloc); */

	/* due to 1-based indexing in LUmod we need to shift vectors backwards */
	mode = 1; /* expand */
	LUmod(mode, L->rows_alloc, L->rows, tmp, tmp, L->A-1, U->A-1, y-1, z-1, w-1);

}

/* Replace the kth column of X */
/*         On entry: */
/*               kcol        says which column of C is being replaced. */
/*               z(*)        contains the new column. */
/*         On exit: */
/*               y(*), w(*)  are altered. */
/*               z(*)        is not altered. */
void LU_Replace_Col(PT_Matrix L, PT_Matrix U, ptrdiff_t kcol, double *y, double *z, double *w)
{
	int mode;
	ptrdiff_t tmp=0;
	kcol++; /* Fortran is 1-indexed */

	/* due to 1-based indexing in LUmod we need to shift vectors backwards */
	mode = 2; /* replace column */
	LUmod(mode, L->rows_alloc, L->rows, tmp, kcol, L->A-1, U->A-1, y-1, z-1, w-1);

}

/* Replace the kth row of X */
/*         On entry: */
/*               krow        says which row of C is being replaced. */
/*               y(*)        contains the new row. */
/*         On exit: */
/*               y(*), z(*), w(*)  are altered. */
void LU_Replace_Row(PT_Matrix L, PT_Matrix U, ptrdiff_t krow, double *y, double *z, double *w)
{
	int mode;
	ptrdiff_t tmp=0;
	krow++; /* Fortran is 1-indexed */

	/* due to 1-based indexing in LUmod we need to shift vectors backwards */
	mode = 3; /* replace row */
	LUmod(mode, L->rows_alloc, L->rows, krow, tmp, L->A-1, U->A-1, y-1, z-1, w-1);

}

/* Delete the krow-th row and kcol-th column */
/*         On entry: */
/*               krow, kcol  say which row and col of C are being replaced. */
/*         On exit: */
/*               y(*), z(*), w(*)  are altered. */
void LU_Shrink(PT_Matrix L, PT_Matrix U, ptrdiff_t krow, ptrdiff_t kcol, double *y, double *z, double *w)
{
	int mode;
	krow++; /* Fortran is 1-indexed */
	kcol++;

	/* due to 1-based indexing in LUmod we need to shift vectors backwards */
	mode = 4; /* shrink */
	LUmod(mode, L->rows_alloc, L->rows, krow, kcol, L->A-1, U->A-1, y-1, z-1, w-1);
	L->rows = L->rows - 1;
	L->cols = L->cols - 1;
	U->rows = U->rows - 1;
	U->cols = U->cols - 1;

}

/************************************
  Given the factorization LB = U for some B, solve the problem
  Bx = vec for x
  Solve using LUMOD functions.
************************************/
void LU_Solve0(PT_Matrix pL, PT_Matrix pU, double *vec, double *x)
{
	int mode;
	ptrdiff_t n;

	n = Matrix_Rows(pL);

	/* solve using lumod */
	/* solve for Bx = vec */
	mode = 1;

	/* due to 1-based indexing in Lprod, Usolve we need to shift vectors backwards */
	/* Computes x = L*vec */
	Lprod(mode, pL->rows_alloc, n, pL->A-1, vec-1, x-1);

	/* Computes x_new s.t. U x_new = x */
	Usolve(mode, pU->rows_alloc, n, pU->A-1, x-1);
	
	/* Vector_Print_raw(x,n); */
}


/************************************
  Given the factorization LB = U for some B, solve the problem
  Bx = vec for x
  Solve using LAPACK functions.
************************************/
void LU_Solve1(PT_Matrix pL, PT_Matrix pUt, double *vec, double *x, ptrdiff_t *info)
{
	ptrdiff_t n, incx=1;
	char U='U', N='N',T='T';
	double alpha=1.0, beta=0.0;

	n = Matrix_Rows(pL);

	/* solve using lapack */
	/* compute x = L*vec */
	dgemv(&T, &n, &n, &alpha, pL->A, &(pL->rows_alloc), vec, &incx, &beta, x, &incx);

	/* solve U*xnew = x using lapack function that also checks for singularity */
	dtrtrs(&U, &N, &N, &n, &incx, pUt->A, &(pUt->rows_alloc), x, &n, info);

	/* printf("my x=\n"); */
	/* Vector_Print_raw(x,n); */
}


void LU_Refactorize(PT_Basis pB)
{
	char L = 'L'; /* lower triangular */
	char D = 'U'; /* unit triangular matrix (diagonals are ones) */
	ptrdiff_t info, incx=1, incp;
	
	/* Matrix_Print_row(pB->pLX); */
	/* Matrix_Print_utril_row(pB->pUX); */

	/* factorize using lapack */
	dgetrf(&(Matrix_Rows(pB->pF)), &(Matrix_Rows(pB->pF)),
	       pMAT(pB->pF), &((pB->pF)->rows_alloc), pB->p, &info);

	/* store upper triangular matrix (including the diagonal to Ut), i.e. copy Ut <- F */
	/* lapack ignores remaining elements below diagonal when computing triangular solutions */
	Matrix_Copy(pB->pF, pB->pUt, pB->w);

	/* transform upper part of F (i.e. Ut) to triangular row major matrix UX*/
	/* UX <- F */
	Matrix_Full2RowTriangular(pB->pF, pB->pUX, pB->r);

	/* invert lower triangular part  */
	dtrtri( &L, &D, &(Matrix_Rows(pB->pF)), pMAT(pB->pF),
		&((pB->pF)->rows_alloc), &info);
			
	/* set strictly upper triangular parts to zeros because L is a full matrix
	 * and we need zeros to compute proper permutation inv(L)*P */
	Matrix_Uzeros(pB->pF);

	/* transpose matrix because dlaswp operates rowwise  and we need columnwise */
	/* LX <- F' */
	Matrix_Transpose(pB->pF, pB->pLX, pB->r);

	/* interchange columns according to pivots in pB->p and write to LX*/
	incp = -1; /* go backwards */
	dlaswp( &(Matrix_Rows(pB->pLX)), pMAT(pB->pLX), &((pB->pLX)->rows_alloc),
		&incx, &(Matrix_Rows(pB->pLX)) , pB->p, &incp);

	/* Matrix_Print_col(pB->pX); */
	/* Matrix_Print_row(pB->pLX); */
	/* Matrix_Print_col(pB->pUt); */
	/* Matrix_Print_utril_row(pB->pUX); */

	/* matrix F after solution is factored in [L\U], we want the original format for the next call
	   to dgesv, thus create a copy F <- X */
	Matrix_Copy(pB->pX, pB->pF, pB->w);


}

/************************************************************************
  Matrix functions 
************************************************************************/

PT_Matrix Matrix_Init(ptrdiff_t m, ptrdiff_t n, const char *name)
{
	PT_Matrix pMat;
	pMat       = (PT_Matrix)malloc(sizeof(T_Matrix));
	pMat->A    = (double*)malloc(m*n*sizeof(double));
	pMat->name = (char*)malloc((strlen(name)+1)*sizeof(char));
	strcpy(pMat->name,name);
	pMat->rows = m;
	pMat->cols = n;
	pMat->rows_alloc = m;
	pMat->cols_alloc = n;

	return pMat;
}

PT_Matrix Matrix_Init_NoAlloc(ptrdiff_t m, ptrdiff_t n, const char *name)
{
	PT_Matrix pMat;
	pMat       = (PT_Matrix)malloc(sizeof(T_Matrix));
	pMat->name = (char*)malloc((strlen(name)+1)*sizeof(char));
	strcpy(pMat->name,name);
	pMat->rows = m;
	pMat->cols = n;
	pMat->rows_alloc = m;
	pMat->cols_alloc = n;

	return pMat;
}


void Matrix_Free(PT_Matrix pMat)
{
	free(pMat->A);
	free(pMat->name);
}

void Matrix_Free_NoAlloc(PT_Matrix pMat)
{
	free(pMat->name);
}

void print_num(double x)
{
	printf("%9.4f ",x);
}


/* Adds the row to the last rows of pA */
/* Assumes column-major format */
void Matrix_AddRow(PT_Matrix pA, double *row)
{
	ptrdiff_t incy = 1;
  
/*  ASSERT(pA->rows < pA->rows_alloc)*/
  
	dcopy(&(Matrix_Cols(pA)),
	       row, &incy,
	       &(C_SEL(pA,pA->rows,0)), &(pA->rows_alloc));
	pA->rows = pA->rows + 1;
}

/* Appends the column to pA */
/* Assumes column-major */
void Matrix_AddCol(PT_Matrix pA, double *col)
{
	ptrdiff_t incx = 1;
  
/*  ASSERT(pA->cols < pA->cols_alloc)*/

	dcopy(&(Matrix_Rows(pA)),
	       col, &incx,
	       &(C_SEL(pA,0,pA->cols)), &incx);
	pA->cols = pA->cols + 1;
}

/* Moves last row of pA to row and then removes last row */
/* Assumes column-major */
void Matrix_RemoveRow(PT_Matrix pA, ptrdiff_t row)
{
/*  ASSERT(pA->rows > 0)*/

	dcopy(&(Matrix_Cols(pA)),
	       &(C_SEL(pA,Matrix_Rows(pA)-1,0)), &(pA->rows_alloc),
	       &(C_SEL(pA,row,0)), &(pA->rows_alloc));
	pA->rows = pA->rows - 1;
}

/* Moves last column of pA to col and then removes last row */
/* Assumes column-major */
void Matrix_RemoveCol(PT_Matrix pA, ptrdiff_t col)
{
	ptrdiff_t incx = 1;
/*  ASSERT(pA->cols > 0)*/

	dcopy(&(Matrix_Rows(pA)),
	       &(C_SEL(pA,0,Matrix_Cols(pA)-1)), &incx,
	       &(C_SEL(pA,0,col)), &incx);
	pA->cols = pA->cols - 1;
}

/* Replaces a column with the given vector */
void Matrix_ReplaceCol(PT_Matrix pA, ptrdiff_t kcol, double *vec)
{
	ptrdiff_t incx = 1;
/*  ASSERT(kcol < pA->cols)*/

	dcopy(&(Matrix_Rows(pA)),
	       vec, &incx,
	       &(C_SEL(pA,0,kcol)), &incx);
}

/* Replaces a row with the given vector */
void Matrix_ReplaceRow(PT_Matrix pA, ptrdiff_t krow, double *vec)
{
	ptrdiff_t incx = 1;
/*  ASSERT(krow < pA->rows)*/

	dcopy(&(Matrix_Cols(pA)),
	       vec, &incx,
	       &(C_SEL(pA,krow,0)), &(pA->rows_alloc));
}

/* Copies matrix B <- A column-wise */
void Matrix_Copy(PT_Matrix pA, PT_Matrix pB, double *w)
{
	ptrdiff_t i, j;

	memcpy(&(pB->rows), &(Matrix_Rows(pA)), 1*sizeof(ptrdiff_t));
	memcpy(&(pB->cols), &(Matrix_Cols(pA)), 1*sizeof(ptrdiff_t));
	memcpy(&(pB->rows_alloc), &(pA->rows_alloc), 1*sizeof(ptrdiff_t));
	memcpy(&(pB->cols_alloc), &(pA->cols_alloc), 1*sizeof(ptrdiff_t));

	for (i=0;i<Matrix_Cols(pA);i++) {
		for (j=0;j<Matrix_Rows(pA);j++) {
			w[j] = C_SEL(pA, j, i);
		}
		Matrix_ReplaceCol(pB, i, w);
	}
}

/* transposes matrix B <- A' column-wise */
void Matrix_Transpose(PT_Matrix pA, PT_Matrix pB, double *w)
{
	ptrdiff_t i, j;

	memcpy(&(pB->cols), &(Matrix_Rows(pA)), 1*sizeof(ptrdiff_t));
	memcpy(&(pB->rows),&(Matrix_Cols(pA)), 1*sizeof(ptrdiff_t));
	memcpy(&(pB->cols_alloc), &(pA->rows_alloc), 1*sizeof(ptrdiff_t));
	memcpy(&(pB->rows_alloc), &(pA->cols_alloc), 1*sizeof(ptrdiff_t));

	for (i=0; i<Matrix_Cols(pA); i++) {
		for (j=0; j<Matrix_Rows(pA); j++) {
			w[j] = C_SEL(pA, j, i);
		}
		Matrix_ReplaceRow(pB, i, w);
	}

}

/* creates column major upper triangular matrix F from U */
void Matrix_Triu(PT_Matrix pF, PT_Matrix pU, double *w)
{
	ptrdiff_t r, c;
	/* creates a copy F <- U */
	Matrix_Copy(pU, pF, w);
	/* put upper diagonal elements to proper position in matrix F */
	for (r=0; r<pU->rows; r++) {
		for (c=0; c<pU->cols; c++) {
			if (c<r) {
				/* lower triangular parts of pF are zeros */
				pF->A[r+c*pF->rows_alloc] = 0.0;		
			}
			else {
				/* access upper triangular elements from matrix pU */
				pF->A[r+c*pF->rows_alloc] = pU->A[pU->rows_alloc*r+c-(r+1)*r/2];
			}
		}
	}
}

/* creates row major triangular UX matrix from a full triangular matrix F */
void Matrix_Full2RowTriangular(PT_Matrix pF, PT_Matrix  pU, double *w)
{
	ptrdiff_t r, c;
	/* transpose  U <- F' */
	Matrix_Transpose(pF, pU, w);
	/* put upper diagonal elements to proper position in matrix F */
	for (r=0; r<pU->rows; r++) {
		for (c=0; c<pU->cols; c++) {
			if (c>=r) {
				/* access upper triangular elements from matrix pU */
				pU->A[pU->rows_alloc*r+c-(r+1)*r/2] = pF->A[r+c*pF->rows_alloc];
			}
		}
	}
}

/* set strictly upper triangular elements as zeros and diagonal as ones*/
void Matrix_Uzeros(PT_Matrix pF)
{
	ptrdiff_t r, c;
	
	for (r=0; r<pF->rows; r++) {
		for (c=0; c<pF->cols; c++) {
			if (c==r) {
				pF->A[r+c*pF->rows_alloc] = 1.0;		
			}
			if (c>r) {
				pF->A[r+c*pF->rows_alloc] = 0.0;		
			}
		}
	}
}


/* Print out a matrix */
/* Column-major formulation */
void Matrix_Print_col(PT_Matrix pMat)
{
	ptrdiff_t r,c;
    
	printf("%s = [\n",pMat->name);
	for(r=0;r<pMat->rows;r++)
	{
		printf("  ");
		for(c=0;c<pMat->cols;c++)
			print_num(C_SEL(pMat,r,c));
		printf("\n");
	}
	printf("];\n");
}

/* Print out a matrix */
/* Row-major formulation */
void Matrix_Print_row(PT_Matrix pMat)
{
	ptrdiff_t r,c;
	double *mat;
	mat = pMat->A;
          
	printf("%s = [\n",pMat->name);
	for(r=0;r<pMat->rows;r++)
	{
		printf("  ");
		for(c=0;c<pMat->cols;c++)
			print_num(R_SEL(pMat,r,c));
		printf("\n");
	}
	printf("];\n");
}

/* Print and upper-triangular matrix */
void Matrix_Print_utril_row(PT_Matrix pMat)
{
	ptrdiff_t r,c,k;
	double *mat;
	mat = pMat->A;
            
	printf("%s = [\n",pMat->name);
	k=0;
	for(r=0;r<pMat->rows;r++)
	{
		printf("  ");
		for(c=0;c<r;c++)
			print_num(0.0);
		for(;c<pMat->cols;c++)
			print_num(mat[k++]);
		for(;c<pMat->cols_alloc;c++) k++;
		printf("\n");
	}
	printf("];\n");
}

void Matrix_Print_size(PT_Matrix pMat)
{
	printf("size(%s) = %ld x %ld, allocated size = %ld x %ld\n", 
	       pMat->name, pMat->rows, pMat->cols,
	       pMat->rows_alloc, pMat->cols_alloc);
}

void Vector_Print_raw(double *vec, ptrdiff_t m)
{
	ptrdiff_t i;
    
	printf("Vector = [\n");
	for(i=0;i<m;i++)
	{
		print_num(vec[i]);
		printf(" ");
	}
	printf("];\n");
}


/* Copy the row pA_rind,pI into row  */
/* pA assumed column major */
void Matrix_GetRow(PT_Matrix pA, double *row, ptrdiff_t row_index, PT_Index col_index)
{
	ptrdiff_t i;
	for(i=0;i<Index_Length(col_index);i++)
		row[i] = C_SEL(pA,row_index,Index_Get(col_index,i));
}

void Matrix_GetCol(PT_Matrix pA, double *col, PT_Index row_index, ptrdiff_t col_index)
{
	ptrdiff_t i;
	for(i=0;i<Index_Length(row_index);i++)
		col[i] = C_SEL(pA, Index_Get(row_index,i), col_index);
}



/************************************************************************
 Index functions
************************************************************************/

PT_Index Index_Init(ptrdiff_t len, const char *name)
{
	PT_Index pIndex;
	pIndex       = (PT_Index)malloc(sizeof(T_Index));
	pIndex->I    = (ptrdiff_t*)malloc(len*sizeof(ptrdiff_t));
	pIndex->name = (char*)malloc((strlen(name)+1)*sizeof(char));
	strcpy(pIndex->name,name);
	pIndex->len  = len;
	pIndex->len_alloc = len;

	return pIndex;
}

void Index_Free(PT_Index pIndex)
{
	free(pIndex->I);
	free(pIndex->name);
}

/* Print out an index */
void Index_Print(PT_Index pI)
{
	ptrdiff_t i;

	printf("%s = [", pI->name);
	for(i=0;i<pI->len;i++)
		printf("%ld ",pI->I[i]);
	printf("]\n");
}

/* Remove the leave'th element (not the element that equals leave) */
void Index_Remove(PT_Index pI, ptrdiff_t leave)
{
	pI->I[leave] = pI->I[pI->len-1];
	pI->len = pI->len - 1;
}

void Index_Add(PT_Index pI, ptrdiff_t num)
{
	/* ASSERT(pI->len < pI->len_alloc)*/
  
	pI->I[Index_Length(pI)] = num;
	pI->len = pI->len + 1;
}

void Index_Replace(PT_Index pI, ptrdiff_t index, ptrdiff_t newval)
{
	/* ASSERT(pI->len > index); */

	pI->I[index] = newval;
}

/* t = find(pI == element) */
ptrdiff_t Index_FindElement(PT_Index pI, ptrdiff_t element)
{
	ptrdiff_t i;
	for(i=0;i<Index_Length(pI);i++)
		if(Index_Get(pI,i) == element)
			return i;
	return -1;
}
