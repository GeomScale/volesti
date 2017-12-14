/*
  FILE DETAILS
  description: Implementation of lexicographic Lemke's algorithm. 
  literature:
  1) Katta G. Murty, Vincent F. Yu; Linear Complementarity, Linear and Nonlinear Programming, internet
  edition, 1997
  available at: http://ioe.engin.umich.edu/people/fac/books/murty/linear_complementarity_webbook/

  2) Richard W. Cottle, Jong-Shi Pang, Richard E. Stone, The linear complementarity problem, Computer
  science and scientific computing, Academic Press 1992

  project: MPT 3.0
  filename: lcp_main.c
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
  details: Added options structure to control the verbosity and maximum number of pivots. Options
  structure enters Basis_Solve, LexMinRatio functions where other details can passed from the user.
  Added refactorization of basis to avoid numerical problems in lex-min ratio.
  If the lex-pivot fails, the pivot element is taken as maximum from v_i over a set of ties in q_i/v_i.
  Other robust suggestions can be found here: Robust implementation of Lemke's method for the linear
  complementarity problem, J. A. Tomlin, Mathematical Programming Studies, 1978, Volume 7, 55-60,
  DOI: 10.1007/BFb0120781   
*/


/************************************************************************
 MODULES 
************************************************************************/

/* standard C headers are included in "lcp_matrix.h" which is a part of "lcp.h" thus
   only main LCP header is provided */
#include "lcp.h"


/************************************************************************
 FUNCTIONS 
************************************************************************/

/* Prepare the matrices for the lcp call

   pA = [-M -1 q]

*/
PT_Matrix lcp_Matrix_Init(ptrdiff_t m, double *M, double *q)
{
	ptrdiff_t i,n,incx;
	double tmp;

	PT_Matrix pA;
      
	/* Build column-major matrix A = [-M -1 q] */
	pA = Matrix_Init(m,m+2,"A");

	/* A(:,1:m) = M */
	memcpy(pMAT(pA), M, m*m*sizeof(double));

	/* A(:,1:m) = -A(:,1:m) */
	tmp  = -1.0;
	incx = 1;
	n = m*m;
	dscal(&n, &tmp, pMAT(pA), &incx);

	/* A(:,m+1) = -1 */
	for(i=0;i<m;i++)
		C_SEL(pA,i,m) = -1.0;

	/* A(:,m+2) = q */
	memcpy(&(C_SEL(pA,0,m+1)),q,m*sizeof(double));

	return pA;
}


/* Compute the smallest row vector */
/* leave = argmin{iB_i*[q I]/v_i | v_i > 0} */
/* leave is an index into the current basis [W Z+m]*/
/* On exit v is modified */
/* If something goes wrong, returned value contains LCP_flag (negative values). */
ptrdiff_t LexMinRatio(PT_Basis pB, PT_Matrix pA, double *v, T_Options options)
{
	ptrdiff_t m, i, j, k, incx, info;
	PT_Index pP;
	double *rat, *q;
	double minrat=0.0, alpha;

	m   = Matrix_Rows(pA);
	pP  = pB->pP;
	rat = pB->y;
	q   = pB->w;

	/* used to obtain an identity matrix in lex-pivot operation */
	alpha = 0.0;
	/* increment of each vectors is 1 */
	incx = 1;
  
	/* Compute the index set of positive v s */
	Index_Length(pP) = 0;
	for (j=0, i=0;i<m;i++)
	{
		if ( v[i] > options.zerotol )
		{
			Index_Add(pP, i);
			/* The positive elements of 1/v are now the first |pP| elements */ 
			v[j++] = v[i];
		}
	}

	if(Index_Length(pP) == 0) /* Problem infeasible */
		return LCP_INFEASIBLE;

	/* q = A(:,m+1) */
	memcpy(q,&(C_SEL(pA, 0, m+1)),m*sizeof(double));
	
        for (i=0; i<m+1; i++)
	 {	
		/* Basis_Print(pB); */
		/* Index_Print(pP); */

		/* rat = iB*q */
		info = Basis_Solve(pB, q, rat, options);
		if (info!=0)  {
#ifdef MATLAB_MEX_FILE
		  if (options.verbose) {
		    printf("Numerical problems with Basis_Solve in LexRatio.\n");
		  }
#endif
			return LCP_CODEERROR;
		}	

		/* printf("q after Basis_Solve in LexRatio=\n"); */
		/* Vector_Print_raw(q,m); */

		/* rat = rat(pP) / v(pP) */
		/* find minimum ratio */
		minrat = rat[Index_Get(pP,0)]/v[0];
		for (j=0; j<Index_Length(pP); j++) {
			rat[j] = rat[Index_Get(pP,j)]/v[j];
			if (rat[j] < minrat) {
				minrat = rat[j];
			}
		}
		/* printf("minrat=%lf\n",minrat); */
		/* printf("Ratio vector after=\n"); */
		/* Vector_Print_raw(rat,m); */

		/* printf("positive index set:\n"); */
                /* Index_Print(pP); */
    
		/* remove all indices that are greater than the minumum ratio with some tolerance */
		/* Move all elements equal to minrat to the start of the list */
		for (k=0, j=0;j<Index_Length(pP);j++)
		{
			/* printf("j=%ld, fabs(minrat - rat[lndex_Get(pP,j)]) = %f\n", j, fabs(minrat - rat[Index_Get(pP,j)])); */
			if (fabs(minrat - rat[j]) < options.lextol) {
				Index_Set(pP,k,Index_Get(pP,j));
				/*pP is now a set of indices with equal ratio */
				v[k++] = v[j];
			} 
		}
		pP->len = k;
		/* printf("new index set:\n"); */
                /* Index_Print(pP); */

		if(Index_Length(pP) <= 0) {
			/* Vector_Print_raw(v,m); */
			/* Vector_Print_raw(rat,m); */
			/* mexPrintf("minrat = %lf, length(pP)  = %ld\n",minrat,Index_Length(pP)); */
			/* Index_Print(pP); */
            #ifdef MATLAB_MEX_FILE
		  if (options.verbose) {
			printf("Numerical problems with finding lexicographic minimum. You can refine this using the 'lextol' feature in options.\n");
		  }
            #endif
			return LCP_CODEERROR;
		}

		if(Index_Length(pP) > 1)
		{
			/* printf("Taking lex-pivot\n"); */
			/* Need to make a lex-comparison */
			/* q = 0*q */
			dscal(&m, &alpha, q, &incx);
			q[i] = 1.0;
		}
		else
		{
			return Index_Get(pP,0);
		}
	 }

	/* if the lex-pivot failed, pick the largest element from v_i vector i=0, ..., k */
	return Index_Get(pP, idamax(&k, v, &incx)-1);

	/* mexErrMsgTxt("Did not find a lex-min. This should be impossible."); */
	/* return LCP_INFEASIBLE; */
}


  
/*
  The main LCP function

  Goal: compute w,z > 0 s.t. w-Mz = q, w'z = 0
  
  On Input:

  pA = [-M -1 q]

  Returns:
  LCP_FEASIBLE    1
  LCP_INFEASIBLE -1
  LCP_UNBOUNDED  -2
  LCP_PRETERMINATED -3
  LCP_CODEERROR  -4

*/
int lcp(PT_Basis pB, PT_Matrix pA, ptrdiff_t *pivs, T_Options options)
{
	ptrdiff_t enter, leave, left;
	ptrdiff_t m, i, info;
	int FirstLoop;
	double *v, *t, total_time;

#ifdef  MATLAB_MEX_FILE
	clock_t t1,t2;
#endif
  
  
	m = Matrix_Rows(pA);
	v = pB->v;
	t = pB->y;
  
	/* Index of the artificial variable */
	enter = 2*m;

#ifdef  MATLAB_MEX_FILE
	t1 = clock();
#endif
	FirstLoop = 1;
	(*pivs) = 0;
	while(1)
	{
		if ( ++(*pivs) >= options.maxpiv ) {
            #ifdef MATLAB_MEX_FILE
		  if (options.verbose) {
			printf("Exceeded maximum number of pivots, returning current basis.\n");
		  }
            #endif
			/* take out the artificial variable from the basis */
			left = Basis_Pivot(pB, pA, enter, Index_Length(pB->pW));
			return LCP_PRETERMINATED;
		}

        #ifdef  MATLAB_MEX_FILE
		t2 = clock();
		total_time = ((double)(t2-t1))/CLOCKS_PER_SEC;
		#else
		total_time = -1;
		#endif

		if ( total_time >= options.timelimit ) {
            #ifdef MATLAB_MEX_FILE
		  if (options.verbose) {
			printf("Maximum time limit was achieved, returning current basis.\n");
		  }
            #endif
			/* take out the artificial variable from the basis */
			left = Basis_Pivot(pB, pA, enter, Index_Length(pB->pW));
			return LCP_PRETERMINATED;
		}


		/******************************/
		/* 1. Choose leaving variable */
		/******************************/

		if(FirstLoop == 1)
		{
			FirstLoop = 0;
			for(i=0;i<m;i++) v[i]=1.0;
		} else
		{
			/* Solve for v = inv(B)*Ae */
			if(enter < m)
			{
				/* A(:,enter) is [0 ... 1 ... 0] */
				for(i=0;i<m;i++) t[i]=0.0;
				t[enter] = 1.0;
			}
			else
			{
				/* v = M(:,enter-m) */
				memcpy(t,&(C_SEL(pA,0,enter-m)),sizeof(double)*m);
			}

            #ifdef MATLAB_MEX_FILE
			if (options.verbose) {
				printf("right hand side entering basis_solve:\n");
				Vector_Print_raw(t,m);
			}
            #endif

			info = Basis_Solve(pB, t, v, options);
			if (info!=0)  {
                #ifdef MATLAB_MEX_FILE
			  if (options.verbose) {
				printf("info=%ld\n", info);
				printf("Numerical problems in lcp.\n");
			  }
                #endif
				return LCP_INFEASIBLE;
			}	


            #ifdef MATLAB_MEX_FILE
			/* print basis if required */
			if (options.verbose) {
				Basis_Print(pB);
				printf("Basic solution:\n");
				Vector_Print_raw(v,m);
			}
            #endif

		}

		/* printf("vector v:\n"); */
		/* Vector_Print_raw(v, m); */
		leave = LexMinRatio(pB, pA, v, options);
		/* Test if leave <0  => no valid leaving var, returning LCP_flag */
		if (leave < 0) {
			return leave;
		}
        
        #ifdef MATLAB_MEX_FILE
		if (options.verbose) {
			printf("Lexicographic minimum found with leaving var=%ld\n",leave); 
		}
        #endif
   
    
		/******************************/
		/* 2. Make the pivot          */
		/******************************/

        #ifdef MATLAB_MEX_FILE
		if (options.verbose) {
			printf("Pivoting with variables: enter = %4ld leave = %4ld ...\n", enter, leave);
		}
        #endif

		/* do the pivot */
		left = Basis_Pivot(pB, pA, enter, leave);

		/* for LUMOD we need to refactorize at every n-steps due numerical problems */
		/* every n-steps refactorize L, U from X using lapack's routine dgetrf */
		if  ((options.routine==0) && (*(pivs) % options.nstepf == 0) ) {
			LU_Refactorize(pB);
		}
        
        #ifdef MATLAB_MEX_FILE
		if (options.verbose) {
			printf("Leaving variable after pivot operation: left = %4ld\n", left);
		}
        #endif
    
		/* Did the artificial variable leave? */
		if(left == 2*m) {
			return LCP_FEASIBLE;
		}

		/*******************************/
		/* 3. Choose entering variable */
		/*******************************/
    
		/* The entering variable must be the
		   complement of the previous leaving variable */
		enter = (left+m) % (2*m);
	}
}
