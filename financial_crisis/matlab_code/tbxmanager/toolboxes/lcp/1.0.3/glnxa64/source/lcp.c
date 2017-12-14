/* 
FILE DETAILS
description: Matlab interface to lexicographic Lemke algorithm.
project: MPT 3.0
filename: lcp.c


COMPILATION:
Requires BLAS, LAPACK libraries installed (should be by default in Matlab).
Requires lcp.h, lcp_main.c, lcp_matrix.h, lcp_matrix.c, lumod_dense.h, lumod_dense.c
To compile this to Matlab executable file under LINUX/MAC, use the following syntax:

    mex -largeArrayDims lcp.c lcp_main.c lcp_matrix.c lumod_dense.c  -lmwblas -lmwlapack

or call "lcp_compile" from Matlab.

 PURPOSE:
  Matlab interface to lexicographic Lemke algorithm which solves linear-complementarity problem
  (LCP) of the form

            find w, z
   s.t.:
          w - M*z = q
               w >= 0
               z >= 0
             w'*z = 0

  
  SYNTAX:
              [z, w, basis, exitflag, pivots, time] = lcp(M, q, options);

  inputs:
          M - positive semi-definite square patrix (for a feasible solution)
          q - right hand side vector with dimension equal M
  
  options:
         .zerotol      (1e-10)  = Less than this treshold the value is considered as zero
         .lextol       (1e-10) = Lexicographic tolerance - a small treshold from which values
                               are considered as equal.
         .maxpiv      (1e6)  = Maximum number of pivots to be performed
         .nstepf        (50) = If routine is 0, then every 50 pivot steps the basis is
                      refactorized to avoid numerical problems for LUMOD. For
                      other routines the factorization is performed at
                      each step.
         .clock       ([0] or 1) = Show the information about the computational time.
         .verbose      ([0] or 1) = Verbose output. Show progress of pivoting algorithm including
                       entering, leaving variables, actual basis and basis solution.
         .routine      ([0] or 1, 2) = Routine which should be used to obtain a basis solution.
                          0 - corresponds to LUmod package that performs
                          factorization in the form L*A = U. Depending on
                          the change in A factors L, U are updated. This is
                          the fastest method. 
                          1 - corresponds to DGESV simple driver from LAPACK
                          package which solves the system AX = B by
                          factorizing A and overwriting B with the solution
                          X. Since the factorization is performed at each
                          pivot step, this method tends to be much slower
                          than method 0.
                          2 - corresponds to DGELS simple driver 
                          which solves overdetermined or underdetermined
                          real linear systems min ||b - Ax||_2 involving an
                          M-by-N matrix A, or its transpose, using a QR or
                          LQ  factorization of A. Since the factorization
                          is performed at each pivot step, this method tends
                          to be much slower than method 0.
         .timelimit    (3600)  = Time limit in seconds. If this limit is exceeded, the pivoting
                              algorithm is terminated and current basis is returned.
         .normalize (0 or [1]) = Input matrices M, q get scaled by D1 and D2 when invoking this
                              option: Mn = D1*M*D2, qn = D1*q, and solution is recovered as 
                              z = D2*zn, w = M*z+q
         .normalizethres (1e6) = If the normalize option is on, then the matrix scaling is performed 
                              only if 1 norm of matrix M (maximum absolute column sum) is above this
                              threshold. 

 
 outputs:
         z - solution vector to the LCP problem
         w - complementary vector to z
         basis - an index set describing the feasible basis
         exitflag = 1 - feasible solution
                    -1 - infeasible
                    -2 - unbounded
                    -3 - preterminated (due to time limit or maximum pivot limit )
                    -4 - other (numerical) error
         pivots - total number of pivots performed by the algorithm
         time - time needed for algorithm to terminate (in seconds)


AUTHOR:
 Copyright (C) 2006 by Colin N. Jones, Automatic Control Laboratory, ETH Zurich,
 cjones@control.ee.ethz.ch
 
 Revised 2010 by Martin Herceg, Automatic Control Laboratory, ETH Zurich,  
 herceg@control.ee.ethz.ch  


REVISION HISTORY:
  date: August 2011
  details: Corrected detection of the initial feasibility, instead of qn[i] < zerotol 
          it is checked qn[i] < -zerotol
 
  date: August 2011
  details: Added "normalizethres" option which represents a threshold from which to 
         to perform normalization. It is valid only if "normalize" option is on.
 
  date: April 2011
  details: Modified "NormalizeMatrix" function to compute properly row/column norms.
  
  date: April 2011 
  details: Testing if the problem is not feasible at the beginning.
 
  date: March 2011
  details: Minor comments.

  date: Dec, 2010
  details: 
  Default options can be modified by the user.
  Added checks for input arguments.
  Added scaling function which is controllable by options.normalize.


LICENSE:
  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  
*/

#include "lcp.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	ptrdiff_t i, j, m, pivs, inc=1, info;
	int Result, nfields=0, qinfeas=0;
	const char *fname;
	char msg[100], T='n';
	double *w, *z, *I, *x, *r, *c, *Mn, *qn;
	mxArray *fval;
	double total_time, alpha=1.0, s=0.0, sn=0.0;
	clock_t t1,t2;
	T_Options options;  /* options structure defined in lcp_matrix.h */    
	PT_Matrix pA; /* Problem data A = [M -1 q] */
	PT_Basis  pB; /* The basis */  

	/* inputs:
	 * M - positive semi-definite square patrix (for a feasible solution)
	 * q - right hand side vector with dimension equal M
	 * options - structure with options
	 */

	/* tests for input arguments */
	if ( nrhs < 2)
		mexErrMsgTxt("lcp: At least 2 arguments are required.");
	if ( nrhs > 3)
		mexErrMsgTxt("lcp: Too many input arguments.");
	if ( nlhs > 6)
		mexErrMsgTxt("lcp: Too many output arguments.");
    if (mxIsEmpty(prhs[0]))
        mexErrMsgTxt("lcp: First argument must not be empty.");
    if (mxIsEmpty(prhs[1]))
        mexErrMsgTxt("lcp: Second argument must not be empty.");    
	if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) )
		mexErrMsgTxt("lcp: First 2 input arguments must be of type DOUBLE.");
	if ( mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]) )
		mexErrMsgTxt("lcp: Input arguments are not allowed in sparse format.");
	/* only matrices, no tenzors */
	if ( mxGetNumberOfDimensions(prhs[0])!=2 ||
	     mxGetNumberOfDimensions(prhs[1])!=2 )
		mexErrMsgTxt("lcp: First 2 arguments must be matrices."); 
	/* testing options */
	if ( nrhs == 3 ) {
		if (!mxIsStruct(prhs[2]))
			mexErrMsgTxt("lcp: Options must be in STRUCT format.");
		else if (mxGetNumberOfElements(prhs[2])>1)
			mexErrMsgTxt("lcp: Only one option structure is allowed.");
		/* checking each option individually */
		nfields = mxGetNumberOfFields(prhs[2]);
		for (i=0; i<nfields; i++) {
			fname = mxGetFieldNameByNumber(prhs[2], i);
			fval = mxGetField(prhs[2], 0, fname);
			/* check for proper field names */
			if (!( (strcmp(fname, "zerotol")==0) || (strcmp(fname, "maxpiv")==0) ||
			       (strcmp(fname, "lextol")==0) || (strcmp(fname, "nstepf")==0) ||
			       (strcmp(fname, "clock")==0) || (strcmp(fname, "verbose")==0) ||
			       (strcmp(fname, "routine")==0) || (strcmp(fname, "timelimit")==0) ||
			       (strcmp(fname, "normalize")==0) || (strcmp(fname, "normalizethres")==0) )) {
				strcpy(msg,"");
				strcat(msg, "lcp: The field '");
				strcat(msg, fname);
				strcat(msg, "' is not allowed in the options structure.");
				mexErrMsgTxt(msg);
			}
			/* some options must be nonnegative */
			if ( (strcmp(fname,"zerotol")==0) || (strcmp(fname,"maxpiv")==0) ||
			     (strcmp(fname,"timelimit")==0) || (strcmp(fname,"lextol")==0) ||
			     (strcmp(fname,"nstepf")==0) || (strcmp(fname, "normalizethres")==0) ) {
				if (!mxIsDouble(fval) || mxIsEmpty(fval) || (mxGetM(fval)*mxGetN(fval))!=1 ||
				    (mxGetScalar(fval)<=0) ) {	  
					strcpy(msg,"");
					strcat(msg, "lcp: Option value '");
					strcat(msg, fname);
					strcat(msg, "' must be of type double, nonempty, scalar, and nonnegative.");
					mexErrMsgTxt(msg);
				}
			}
			/* some can be zeros */
			if (strcmp(fname,"clock")==0 || strcmp(fname,"verbose")==0 ||
			    strcmp(fname,"routine")==0 || strcmp(fname,"normalize")==0) {
				if (!mxIsDouble(fval) || mxIsEmpty(fval) || (mxGetM(fval)*mxGetN(fval))!=1 ) {	  
					strcpy(msg,"");
					strcat(msg, "lcp: Option value '");
					strcat(msg, fname);
					strcat(msg, "' must be of type double, nonempty, scalar, and integer valued.");
					mexErrMsgTxt(msg);
				}
			}
			/* No NaN values allowed */
			if ( mxIsNaN(*mxGetPr(fval)) )
				mexErrMsgTxt("lcp: No 'NaN' are allowed in the options structure.");
			/* Inf value allowed only for timelimit */
			if ( strcmp(fname,"timelimit")!=0 && mxIsInf(*mxGetPr(fval)) )
				mexErrMsgTxt("lcp: No 'Inf' terms allowed except for 'timelimit' option.");
		}
	}
	/* dimension checks */
	if (mxGetM(prhs[0])!=mxGetN(prhs[0]))
		mexErrMsgTxt("lcp: First input argument must be a square matrix.");
	/* problem dimension */
	m = mxGetM(prhs[0]);
	if ( mxGetNumberOfElements(prhs[1])!=m )	  
		mexErrMsgTxt("lcp: Second input argument must have the same dimension as the first argument.");
	/* No NaN or Inf values allowed */
	for (i=0; i<2; i++) {
		x = mxGetPr(prhs[i]);
		for (j=0; j<mxGetNumberOfElements(prhs[i]); j++) {
			if ( mxIsNaN(x[j]) || mxIsInf(x[j]) )
				mexErrMsgTxt("lcp: No 'NaN' or 'Inf' terms are allowed.");
		}
	}
	
 
	/* default options */
	options.zerotol = 1e-10; /* zero tolerance */
	options.lextol = 1e-10; /* lexicographic tolerance - a small treshold to determine if values are equal */
	options.maxpiv = 1e6; /* maximum number of pivots */
	/* if LUMOD routine is chosen, this options refactorizes the basis after n steps using DGETRF
	   routine from lapack to avoid numerical problems */
	options.nstepf = 50;
	options.clock = 0; /* 0 or 1 - to print computational time */
	options.verbose = 0; /* 0 or 1 - verbose output */
/*  which routine in Basis_solve should solve a set of linear equations: 
    0 - corresponds to LUmod package that performs factorization in the form L*A = U. Depending on
    the change in A factors L, U are updated. This is the fastest method.
    1 - corresponds to DGESV simple driver 
    which solves the system AX = B by factorizing A and overwriting B with the solution X
    2 - corresponds to DGELS simple driver 
    which solves overdetermined or underdetermined real linear systems min ||b - Ax||_2
    involving an M-by-N matrix A, or its transpose, using a QR or LQ  factorization of A.  */
	options.routine = 0; 
	options.timelimit = 3600; /* time limit in seconds to interrupt iterations */
	options.normalize = 1; /* 0 or 1 - perform scaling of input matrices M, q */
    options.normalizethres = 1e6; 
    /* If the normalize option is on, then the matrix scaling is performed 
      only if 1 norm of matrix M (maximum absolute column sum) is above this threshold.
      This enforce additional control over normalization since it seems to be more
      aggressive also for well-conditioned problems. */
      

	if (nrhs==3) {
		/* overwriting default options by the user */
		for(i=0; i<nfields; i++){
			fname = mxGetFieldNameByNumber(prhs[2], i);   
			fval = mxGetField(prhs[2], 0, fname);
			if ( strcmp(fname,"zerotol")==0 )
				options.zerotol = mxGetScalar(fval);
			if ( strcmp(fname,"lextol")==0 )
				options.lextol = mxGetScalar(fval);
			if ( strcmp(fname,"maxpiv")==0 ) {
				if (mxGetScalar(fval)>=(double)INT_MAX)
					options.maxpiv = INT_MAX;
				else
					options.maxpiv = (int)mxGetScalar(fval);
			}
			if ( strcmp(fname,"nstepf")==0 )
				options.nstepf = (int)mxGetScalar(fval);
			if ( strcmp(fname,"timelimit")==0 )
				options.timelimit = mxGetScalar(fval);
			if ( strcmp(fname,"clock")==0 )
				options.clock = (int)mxGetScalar(fval);
			if ( strcmp(fname,"verbose")==0 )
				options.verbose = (int)mxGetScalar(fval);
			if ( strcmp(fname,"routine")==0 )
				options.routine = (int)mxGetScalar(fval);
			if ( strcmp(fname,"normalize")==0 )
				options.normalize = (int)mxGetScalar(fval);		       
            if ( strcmp(fname, "normalizethres")==0 )
                options.normalizethres = mxGetScalar(fval);
		}
	}

	/* Normalize M, q to avoid numerical problems if possible 
	   Mn = diag(r)*M*diag(c) , qn = diag(r)*q  */
	/* initialize Mn, qn */
	Mn = (double *)mxCalloc(m*m, sizeof(double));
	qn = (double *)mxCalloc(m, sizeof(double));
	/* allocate vectors r, c */
	r = (double *)mxCalloc(m, sizeof(double));
	c = (double *)mxCalloc(m, sizeof(double));
	/* initialize to ones */
	for (i=0; i<m; i++) {
		r[i] = 1.0;
		c[i] = 1.0;
	}
	/* write data to Mn = M */
	memcpy(Mn, mxGetPr(prhs[0]), (m*m)*sizeof(double));
	/* write data to qn = q */
	memcpy(qn, mxGetPr(prhs[1]), m*sizeof(double));
    /* check out the 1-norm of matrix M (maximum column sum) */
    for (i=0; i<m; i++) {
        sn = dasum(&m, &Mn[i*m], &inc);
        if (sn>s) {
            s = sn;
        }
    }
    
	/* scale matrix M, q and write scaling factors to r (rows) and c (columns) */        
	if (options.normalize && s>options.normalizethres) {
        NormalizeMatrix (m, m, Mn, qn, r, c,  options);
    }

	/* Setup the problem */
	/* pA = lcp_Matrix_Init(m, mxGetPr(prhs[0]), mxGetPr(prhs[1])); */
	pA = lcp_Matrix_Init(m, Mn, qn);
	pB = Basis_Init(m);

    /* check if the problem is not feasible at the beginning */
    for (i=0; i<m; i++) {
        if (qn[i]<-options.zerotol) {
            qinfeas = 1;
            break;
        }
    }
    
	/* Solve the LCP */
    if (qinfeas) {
        t1 = clock();
        Result = lcp(pB, pA, &pivs, options);
        t2 = clock();
        total_time = ((double)(t2-t1))/CLOCKS_PER_SEC;
    }
    else {
        pivs = 0;
        total_time = 0;
        Result = LCP_FEASIBLE;
    }

	if (options.clock) {
		printf("Time needed to perform pivoting:\n time= %i  (%lf seconds)\n",
		       t2-t1,total_time);
		printf("Pivots: %ld\n", pivs);
		printf("CLOCKS_PER_SEC = %i\n",CLOCKS_PER_SEC);
	}

	/* outputs */
	plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL); /* z */
	plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL); /* w */
	plhs[2] = mxCreateDoubleMatrix(m,1,mxREAL); /* base */
	plhs[3] = mxCreateDoubleScalar(Result); /* exitflag */
	plhs[4] = mxCreateDoubleScalar(pivs); /* pivots */
	plhs[5] = mxCreateDoubleScalar(total_time); /* time */
  
	/* get pointers to outputs */
	z = mxGetPr(plhs[0]);
	w = mxGetPr(plhs[1]);
	I = mxGetPr(plhs[2]);
	/* initialize values to 0 */
	for(i=0;i<m;i++)
	{
		w[i] = 0.0;
		z[i] = 0.0;
		I[i] = 0.0;
	}

	/* for a feasible basis, compute the solution */
	if ( Result == LCP_FEASIBLE || Result == LCP_PRETERMINATED )
	{
		x = (double *)mxCalloc(m, sizeof(double));
		t1 = clock();
		info = Basis_Solve(pB, &(C_SEL(pA,0,m+1)), x, options);
		for (j=0,i=0;i<Index_Length(pB->pW);i++,j++)
		{
			w[Index_Get(pB->pW,i)] = x[j];
			/* add 1 due to matlab 1-indexing */
			I[j] = Index_Get(pB->pW,i)+1;
		}
		for(i=0;i<Index_Length(pB->pZ);i++,j++)
		{
			/* take only positive values */
			if (x[j] > options.zerotol ) {
				z[Index_Get(pB->pZ,i)] = x[j];
			}
			/* add 1 due to matlab 1-indexing */
			I[j] = Index_Get(pB->pZ,i)+m+1;
		}
		t2 = clock();
		if (options.clock) {
			printf("Time in total needed to solve LCP: %lf seconds\n",
			       total_time + (double)(t2-t1)/(double)CLOCKS_PER_SEC);
		}
		
		if (options.normalize) {
			/* do the backward normalization */
			/* z = diag(c)*zn */
			for (i=0; i<m; i++) {
				z[i] = c[i]*z[i];
			}
			
			/* since the normalization does not compute w properly, we recalculate it from
			 * recovered z */
			/* copy w <- q; */
			dcopy(&m, mxGetPr(prhs[1]), &inc, w, &inc);
			/* compute w = M*z + q */
			dgemv(&T, &m, &m, &alpha, mxGetPr(prhs[0]), &m,
			      z, &inc, &alpha, w, &inc);
			/* if w is less than eps, consider it as zero */
			for (i=0; i<m; i++) {
				if (w[i]<options.zerotol) {
					w[i] = 0.0;
				}
			}
		}

		mxFree(x);
	}
  

	/* free allocated variables */
	Matrix_Free(pA);
	Basis_Free(pB); 	
	free(pA); 
	free(pB);

	mxFree(Mn);	
	mxFree(qn);
	mxFree(r);
	mxFree(c);
	
}

/* normalize input matrix M */
void NormalizeMatrix (ptrdiff_t m, ptrdiff_t n, double *A, double *b, double *r, double *c,  T_Options options)
{
    /* scales matrix A by finding D1= diag(r) and D2=diag(c) in An = D1*A*D2  
       such that infinity norm of each row and column approaches 1 

       find D1, D2  
 
       s.t. 
       An = D1*A*D2 

       Scaling matrix is used in solving A*x=b for badly scaled matrix A as
       follows:

       A*x = b       
       D1*A*x = D1*b           / multiply from left by D1 
       D1*A*D2*inv(D2)*x = D1*b    / insert D2*inv(D2)  
       D1*A*D2*y = D1*b            / substitute y = inv(D2)*x 
       An*y = bn              / substitue An = D1*A*D2, v = D1*b 

       First solve An*y = bn, then obtain x = D2*y 

       Details of the method are in file drRAL2001034.ps.gz at 
       http://www.numerical.rl.ac.uk/reports/reports.html. 

       arguments:
       n - (input), ptrdiff_t (INTEGER) 
           Number of rows for matrix A and dimension of vector r
       m - (input), ptrdiff_t (INTEGER)
           Number of cols for matrix A and dimension of vector c
       A -  (input/output) (pointer to DOUBLE)
           Input matrix is changed at the output with the normalized matrix.
           It is assumed that matrix A does not contain row or column with all elements
           zeros. If this not the case, no scaling is performed.          
       b - (inptu/output) (pointer to DOUBLE)
           Input vector is normalized at the output
       r - (output) (pointer to DOUBLE)
           Vector with scaling factors for each row, D1 = diag(r)
       c - (output) (pointer to DOUBLE)
           Vector with scaling factors for each column, D2 = diag(c)
     */

    ptrdiff_t i, j, kmax, inc=1,k=0;
    int flag=1;
    double row_norm = 1.0, col_norm = 1.0, nrm;
    double *rn, *rs, *cn, *cs; 

    /* prepare vectors */
    rn = (double *)mxCalloc(n,sizeof(double));
    rs = (double *)mxCalloc(n,sizeof(double));
    cn = (double *)mxCalloc(m,sizeof(double));
    cs = (double *)mxCalloc(m,sizeof(double));

    /* check if A contains row or cols with all zeros */
    /* initialize output vectors */
    for (i=0; i<m; i++) {
        nrm = dasum(&n, &A[i], &m);
        if (nrm<options.zerotol) {
            flag = 0;
		}
        r[i] = 1.0;
    }
    for (i=0; i<n; i++) {
        nrm = dasum(&m, &A[i*n], &inc);
        if (nrm<options.zerotol) {
            flag = 0;
        }
        c[i] = 1.0;
    }
  
    

	if (flag)
	{
		while ((row_norm>options.zerotol) && (col_norm>options.zerotol)) 
		{
			/* loop thru rows */
			for (i=0; i<m; i++) {
				/* compute infinity norm of each row */
				/* find index of a maximum value in this row (1-indexed) */
				kmax = idamax(&n, &A[i], &m);
				/* evaluate the norm */				
				nrm = fabs(A[i+(kmax-1)*m]);
				rs[i] = 1.0/sqrt(nrm);
				/* printf("rows: kmax=%ld, nrm=%f\n",kmax,nrm); */
				/* to compute row_norm*/
				rn[i] = 1 - nrm;
			}
			/* loop thru cols */
			for (i=0; i<n; i++) {
				/* compute infinity norm of each column */
				/* find maximum value in this column */
				kmax = idamax(&m, &A[i*n], &inc);
				/* evaluate the norm */
				nrm = fabs(A[i*m+kmax-1]);
				cs[i] = 1.0/sqrt(nrm);
				/* printf("cols: kmax=%ld, nrm=%f\n",kmax,nrm); */
				/* set  the col_norm */
				cn[i] = 1 - nrm;
			}

			/* scale matrix  A and vector b */
			for (i=0; i<m; i++) {
				for (j=0; j<n; j++) {
					A[i+j*m] = rs[i]*A[i+j*m]*cs[j];
				}
				b[i] = rs[i]*b[i];
				/* update output vector r */
				r[i] = r[i]*rs[i];	
			}			

			/* update output vectors c*/
			for (i=0; i<n; i++) {
				c[i] = c[i]*cs[i];
			}
				
		
			/* compute the norm and column norms */
			kmax = idamax(&m, rn, &inc);
			row_norm = fabs(rn[kmax-1]);
			kmax = idamax(&n, cn, &inc);
			col_norm = fabs(cn[kmax-1]);
			
			if (++k>=1000 ) {
				/* limit the number of iterations */
				break;
			}
		}
	}
    
	/* freeing vectors */
	mxFree(rn);
	mxFree(rs);
	mxFree(cn);
	mxFree(cs);
		
}
