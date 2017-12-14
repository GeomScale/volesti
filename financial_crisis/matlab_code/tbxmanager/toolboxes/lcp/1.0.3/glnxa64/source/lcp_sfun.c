/* 
FILE DETAILS
description: Matlab Simulink interface to lexicographic Lemke algorithm.
project: MPT 3.0
filename: lcp_sfun.c

COMPILATION:
Requires BLAS, LAPACK libraries installed (should be by default in Matlab).
Requires lcp.h, lcp_main.c, lcp_matrix.h, lcp_matrix.c, lumod_dense.h, lumod_dense.c
To compile this to Matlab executable file under LINUX/MAC, use the following syntax:

    mex -largeArrayDims lcp_sfun.c lcp_main.c lcp_matrix.c lumod_dense.c  -lmwblas -lmwlapack

or invoke the compilation script "lcp_compile" from Matlab.


 PURPOSE:
  Simulink interface to lexicographic Lemke algorithm which solves linear-complementarity problem
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
 Copyrigh (C)  2011 by Martin Herceg, Automatic Control Laboratory, ETH Zurich,  
 herceg@control.ee.ethz.ch  

REVISION HISTORY:
  date: 30 May 2013
  details: Replaced "mxCalloc" and "mxFree" for "calloc" and "free" in ssGetPWork pointers 
          because Matlab R2013 crashed when freeing vectors.
 
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

#define S_FUNCTION_NAME lcp_sfun
#define S_FUNCTION_LEVEL 2

/* icluding all necessary files  */
#include "simstruc.h"
#include "lcp.h"

/* elementary definitions */
#define NSTATES(S)  (int_T)*mxGetPr(ssGetSFcnParam(S,0))  /* first parameter gives the number of states */
#define TS(S) *mxGetPr(ssGetSFcnParam(S,1))   /* sample time, second parameter */
#define NINPUTS(S)   NSTATES(S)*NSTATES(S)+NSTATES(S)  /* number of elements of M + q */
#define NOUTPUTS(S)  6  /* z, w, basis, exitflag, pivots, time */

/*====================*
 * S-function methods *
 *====================*/
 
#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
  /* Function: mdlCheckParameters =============================================
   * Abstract:
   *    Validate our parameters to verify they are okay.
   */
  static void mdlCheckParameters(SimStruct *S)
  {
	  int_T i;
	  static char_T msg[100];

/* for RTW we do not need these variables */
#ifdef MATLAB_MEX_FILE
	const char *fname;
	mxArray *fval;
	int_T nfields;
#endif

	  /* Check 1st parameter: number of states */
      {
	      if ( !mxIsDouble(ssGetSFcnParam(S,0)) || 
		   mxIsEmpty(ssGetSFcnParam(S,0)) ||
		   mxGetNumberOfElements(ssGetSFcnParam(S,0))!=1 ||
		   mxGetScalar(ssGetSFcnParam(S,0))<=0 || 
		   mxIsNaN(*mxGetPr(ssGetSFcnParam(S,0))) || 
		   mxIsInf(*mxGetPr(ssGetSFcnParam(S,0))) ||
		   mxIsSparse(ssGetSFcnParam(S,0)) )  {
		      ssSetErrorStatus(S,"lcp_sfun: Number of states must be of type double, not empty, scalar and finite.");
		      return;
	      }
      }
 
      /* Check 2nd parameter: sample time */
      {
	      if ( !mxIsDouble(ssGetSFcnParam(S,1)) || 
		   mxIsEmpty(ssGetSFcnParam(S,1)) ||
		   mxGetNumberOfElements(ssGetSFcnParam(S,1))!=1 ||
		   mxGetScalar(ssGetSFcnParam(S,1))<=0 || 
		   mxIsNaN(*mxGetPr(ssGetSFcnParam(S,1))) || 
		   mxIsInf(*mxGetPr(ssGetSFcnParam(S,1))) ||
		   mxIsSparse(ssGetSFcnParam(S,1)) )  {
		      ssSetErrorStatus(S,"lcp_sfun: Sample time must be of type double, not empty, scalar and finite.");
		      return;
	      }
      }

#ifdef MATLAB_MEX_FILE
      /* Check 3rd parameter: options */
      if (ssGetSFcnParamsCount(S)==3) 
      {
	      if (!mxIsStruct(ssGetSFcnParam(S,2))) {
		      ssSetErrorStatus(S, "lcp_sfun: Options must be in STRUCT format.");
		      return;
	      }      
	      else if (mxGetNumberOfElements(ssGetSFcnParam(S,2))>1) {
		      ssSetErrorStatus(S,"lcp_sfun: Only one option structure is allowed.");
		      return;
	      }
	      /* checking each option individually */
	      nfields = mxGetNumberOfFields(ssGetSFcnParam(S,2));
	      for (i=0; i<nfields; i++) {
		      fname = mxGetFieldNameByNumber(ssGetSFcnParam(S,2), i);
		      fval = mxGetField(ssGetSFcnParam(S,2), 0, fname);
			/* check for proper field names */
			if (!( (strcmp(fname, "zerotol")==0) || (strcmp(fname, "maxpiv")==0) ||
			       (strcmp(fname, "lextol")==0) || (strcmp(fname, "nstepf")==0) ||
			       (strcmp(fname, "clock")==0) || (strcmp(fname, "verbose")==0) ||
			       (strcmp(fname, "routine")==0) || (strcmp(fname, "timelimit")==0) ||
			       (strcmp(fname, "normalize")==0) || (strcmp(fname, "normalizethres")==0) )) {
				strcpy(msg,"");
				strcat(msg, "lcp_sfun: The field '");
				strcat(msg, fname);
				strcat(msg, "' is not allowed in the options structure.");
				ssSetErrorStatus(S, msg);
				return;
			}
			/* some options must be nonnegative */
			if (strcmp(fname,"zerotol")==0 || strcmp(fname,"maxpiv")==0 ||
			    strcmp(fname,"timelimit")==0 || strcmp(fname,"lextol")==0 ||
			     (strcmp(fname,"nstepf")==0) || (strcmp(fname, "normalizethres")==0) ) {
				if (!mxIsDouble(fval) || mxIsEmpty(fval) || (mxGetM(fval)*mxGetN(fval))!=1 ||
				    (mxGetScalar(fval)<=0) ) {	  
					strcpy(msg,"");
					strcat(msg, "lcp_sfun: Option value '");
					strcat(msg, fname);
					strcat(msg, "' must be of type double, nonempty, scalar, and nonnegative.");
					ssSetErrorStatus(S, msg);
					return;
				}
			}
			/* some can be zeros */
			if (strcmp(fname,"clock")==0 || strcmp(fname,"verbose")==0 ||
			    strcmp(fname,"routine")==0 || strcmp(fname,"normalize")==0) {
				if (!mxIsDouble(fval) || mxIsEmpty(fval) || (mxGetM(fval)*mxGetN(fval))!=1 ) {	  
					strcpy(msg,"");
					strcat(msg, "lcp_sfun: Option value '");
					strcat(msg, fname);
					strcat(msg, "' must be of type double, nonempty, scalar, and integer valued.");
					ssSetErrorStatus(S, msg);
					return;
				}
			}
			/* No NaN values allowed */
			if ( mxIsNaN(*mxGetPr(fval)) ) {
				ssSetErrorStatus(S, "lcp_sfun: No 'NaN' are allowed in the options structure.");
				return;
			}
			/* Inf value allowed only for timelimit */
			if ( strcmp(fname,"timelimit")!=0 && mxIsInf(*mxGetPr(fval)) ) {
				ssSetErrorStatus(S, "lcp_sfun: No 'Inf' terms allowed except for 'timelimit' option.");
				return;
			}
	      }
      }
#endif
    
  }
#endif /* MDL_CHECK_PARAMETERS */


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
	int_T npar = ssGetSFcnParamsCount(S);
	/* options are optional as a third parameter */
	if ( npar==2 || npar==3 ) {
		ssSetNumSFcnParams(S, npar);  /* Number of expected parameters */
	}
	else {
		return;
	}

#ifdef MATLAB_MEX_FILE
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
        mdlCheckParameters(S);
        if (ssGetErrorStatus(S) != NULL) {
            return;
        }
    } else {
        return; /* Parameter mismatch will be reported by Simulink */
    }
#endif

    /* no states because we do not need nothing to remember */
    ssSetNumContStates(S, 0); /* the number of continuous states */
    ssSetNumDiscStates(S, 0 ); /* the number of discrete states */

    /* inputs are M, q*/
    if (!ssSetNumInputPorts(S, 2)) return;
    /* input matrix M to port 0*/
    ssSetInputPortWidth(S, 0, NSTATES(S)*NSTATES(S));     
    ssSetInputPortMatrixDimensions(S, 0, NSTATES(S), NSTATES(S));  
    ssSetInputPortDirectFeedThrough(S, 0, 1); /* direct feedtrough for port 0*/

    /* input vector q to port 1*/
    ssSetInputPortWidth(S, 1, NSTATES(S));     
    ssSetInputPortVectorDimension(S, 1, NSTATES(S));
    ssSetInputPortDirectFeedThrough(S, 1, 1); /* direct feedtrough for port 1*/


    /* outputs are z, w, basis, exitflag, pivots, time */
    if (!ssSetNumOutputPorts(S, 6)) return;
    ssSetOutputPortWidth(S, 0, NSTATES(S)); /* output z to port 0 */
    ssSetOutputPortVectorDimension(S, 0, NSTATES(S));

    ssSetOutputPortWidth(S, 1, NSTATES(S)); /* output w to port 1 */
    ssSetOutputPortVectorDimension(S, 1, NSTATES(S));

    ssSetOutputPortWidth(S, 2, NSTATES(S)); /* basis to port 2 */
    ssSetOutputPortVectorDimension(S, 2, NSTATES(S));

    ssSetOutputPortWidth(S, 3, 1); /* exitflag to port 3 */
    ssSetOutputPortWidth(S, 4, 1); /* pivots to port 4 */
    ssSetOutputPortWidth(S, 5, 1); /* computation time to port 5 */

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 7); /* working vectors Mn, qn, r, c, x, pA, pB */
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* use PWorks to store values */
    ssSetSimStateCompliance(S, HAS_NO_SIM_STATE);

    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE );

}


/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    S-function is comprised of only continuous sample time elements
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
	/* fixed sample time passed from options */
	ssSetSampleTime(S, 0, TS(S));
	ssSetOffsetTime(S, 0, 0.0);
	ssSetModelReferenceSampleTimeDefaultInheritance(S);
}


#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 * Initialize states to zeros.
 */
static void mdlInitializeConditions(SimStruct *S)
{
	int_T  i;
	InputRealPtrsType M, q;
	
	/* checks the input matrices for NaN and Inf */ 

	/* checking M */
	M = ssGetInputPortRealSignalPtrs(S,0);
	for (i=0; i<NSTATES(S)*NSTATES(S); i++) {
		if ( !mxIsFinite(*M[i])  ) {
			ssSetErrorStatus(S, "lcp_sfun: No 'NaN' or 'Inf' terms are allowed in input matrix M.");
			return;
		}
	}
	/* checking q */
	q = ssGetInputPortRealSignalPtrs(S,1);
	for (i=0; i<NSTATES(S); i++) {
		if ( !mxIsFinite(*q[i]) ) {
			ssSetErrorStatus(S, "lcp_sfun: No 'NaN' or 'Inf' terms are allowed in input vector q.");
			return;
		}
	}

}

#define MDL_START  
#if defined(MDL_START)
  /* Function: mdlStart =======================================================
   * Abstract:
   *  allocate memory for work vectors */
static void mdlStart(SimStruct *S)
{
	ssGetPWork(S)[0] = calloc(NSTATES(S)*NSTATES(S), sizeof(double)); /* Mn */
	ssGetPWork(S)[1] = calloc(NSTATES(S), sizeof(double)); /* qn */
	ssGetPWork(S)[2] = calloc(NSTATES(S), sizeof(double)); /* r */
	ssGetPWork(S)[3] = calloc(NSTATES(S), sizeof(double)); /* c */
	ssGetPWork(S)[4] = calloc(NSTATES(S), sizeof(double)); /* x */
	ssGetPWork(S)[5] = Matrix_Init(NSTATES(S),NSTATES(S)+2,"A"); /* pA */
	ssGetPWork(S)[6] = Basis_Init(NSTATES(S));  /* pB */
}
#endif /*  MDL_START */

/* free memory at the end of simulation for work vectors */
static void mdlTerminate(SimStruct *S)
{
    int_T i;
 	for ( i=0; i<5; i++ ) {
		if ( ssGetPWork(S)[i] != NULL ) {
			free( ssGetPWork(S)[i] );
		}
	}
 
	Matrix_Free(ssGetPWork(S)[5]);
	free(ssGetPWork(S)[5]);
    
   	Basis_Free(ssGetPWork(S)[6]);
	free(ssGetPWork(S)[6]);
}

/* Function: mdlOutputs =======================================================
 * do the main optimization routine here
 * no discrete states are considered
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	int_T i, j, Result, qinfeas=0;
	real_T  *z = ssGetOutputPortRealSignal(S,0);
	real_T  *w = ssGetOutputPortRealSignal(S,1);
	real_T  *I = ssGetOutputPortRealSignal(S,2);
	real_T  *exitflag = ssGetOutputPortRealSignal(S,3);
	real_T  *pivots = ssGetOutputPortRealSignal(S,4);
	real_T  *time = ssGetOutputPortRealSignal(S,5);

	InputRealPtrsType M = ssGetInputPortRealSignalPtrs(S,0);
	InputRealPtrsType q = ssGetInputPortRealSignalPtrs(S,1);

	ptrdiff_t pivs, info, m, n, inc=1;
	char_T T='N';
	double total_time,  alpha=1.0, tmp=-1.0, *x, *Mn, *qn, *r, *c, s=0.0, sn=0.0; 
	PT_Matrix pA; /* Problem data A = [M -1 q] */
	PT_Basis  pB; /* The basis */  
	T_Options options;  /* options structure defined in lcp_matrix.h */          

/* for RTW we do not need these variables */
#ifdef MATLAB_MEX_FILE
	const char *fname;
	mxArray *fval;
	int_T nfields;
	clock_t t1,t2;
#endif

    
	UNUSED_ARG(tid); /* not used in single tasking mode */
	
	/* default options */
	options.zerotol = 1e-10; /* zero tolerance */
	options.lextol = 1e-10; /* lexicographic tolerance - a small treshold to determine if values are equal */
	options.maxpiv = INT_MAX; /* maximum number of pivots */
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

#ifdef MATLAB_MEX_FILE
	/* overwriting default options by the user */
	if (ssGetSFcnParamsCount(S)==3) {
		nfields = mxGetNumberOfFields(ssGetSFcnParam(S,2));
		for(i=0; i<nfields; i++){
			fname = mxGetFieldNameByNumber(ssGetSFcnParam(S,2), i);   
			fval = mxGetField(ssGetSFcnParam(S,2), 0, fname);
			if ( strcmp(fname,"zerotol")==0 )
				options.zerotol = mxGetScalar(fval);
			if ( strcmp(fname,"lextol")==0 )
				options.lextol = mxGetScalar(fval);
			if ( strcmp(fname,"maxpiv")==0 ) {
				if (mxGetScalar(fval)>=(double)INT_MAX)
					options.maxpiv = INT_MAX;
				else
					options.maxpiv = (int_T)mxGetScalar(fval);
			}
			if ( strcmp(fname,"nstepf")==0 )
				options.nstepf = (int_T)mxGetScalar(fval);
			if ( strcmp(fname,"timelimit")==0 )
				options.timelimit = mxGetScalar(fval);
			if ( strcmp(fname,"clock")==0 )
				options.clock = (int_T)mxGetScalar(fval);
			if ( strcmp(fname,"verbose")==0 )
				options.verbose = (int_T)mxGetScalar(fval);
			if ( strcmp(fname,"routine")==0 )
				options.routine = (int_T)mxGetScalar(fval);
			if ( strcmp(fname,"normalize")==0 )
				options.normalize = (int_T)mxGetScalar(fval);
            if ( strcmp(fname, "normalizethres")==0 )
                options.normalizethres = mxGetScalar(fval);          
		}
	}
#endif
    

	/* Normalize M, q to avoid numerical problems if possible 
	   Mn = diag(r)*M*diag(c) , qn = diag(r)*q  */
	/* initialize Mn, qn */
	Mn = (double *)ssGetPWork(S)[0];
	qn = (double *)ssGetPWork(S)[1];
	/* initialize vectors r, c */
	r = (double *)ssGetPWork(S)[2];
	c = (double *)ssGetPWork(S)[3];
	/* initialize auxiliary vector x */
	x = (double *)ssGetPWork(S)[4];
	/* initialize to ones */
	for (i=0; i<NSTATES(S); i++) {
		r[i] = 1.0;
		c[i] = 1.0;
	}
	m = NSTATES(S);
	n = m*m;
	/* write data to Mn = M */
	memcpy(Mn, *M, n*sizeof(double));
	/* write data to qn = q */
	memcpy(qn, *q, m*sizeof(double));
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
	pA = ssGetPWork(S)[5];
	/* A(:,1:m) = M */
	memcpy(pMAT(pA), Mn, n*sizeof(double));

	/* A(:,1:m) = -A(:,1:m) */
	dscal(&n, &tmp, pMAT(pA), &inc);

	/* A(:,m+1) = -1 */
	for(i=0;i<m;i++)
		C_SEL(pA,i,m) = -1.0;

	/* A(:,m+2) = q */
	memcpy(&(C_SEL(pA,0,m+1)),qn,m*sizeof(double));

	/* initialize basis */
	pB = ssGetPWork(S)[6];

    /* check if the problem is not feasible at the beginning */
    for (i=0; i<m; i++) {
        if (qn[i]<-options.zerotol) {
            qinfeas = 1;
            break;
        }
    }


    /* Solve the LCP */
    if (qinfeas) {
#ifdef MATLAB_MEX_FILE
	t1 = clock();
#endif

    /* main LCP rouinte */    
    Result = lcp(pB, pA, &pivs, options);
        
#ifdef MATLAB_MEX_FILE
    t2 = clock();
    total_time = ((double)(t2-t1))/CLOCKS_PER_SEC;
#else
    total_time = -1;
#endif
    } else {
        pivs = 0;
        total_time = 0;
        Result = LCP_FEASIBLE;    
    }

#ifdef MATLAB_MEX_FILE
	if (options.clock) {
		printf("Time needed to perform pivoting:\n time= %i  (%lf seconds)\n",
		       t2-t1,total_time);
		printf("Pivots: %ld\n", pivs);
		printf("CLOCKS_PER_SEC = %i\n",CLOCKS_PER_SEC);
	}
#endif

	/* initialize values to 0 */
	for(i=0;i<NSTATES(S);i++)
	{
		w[i] = 0.0;
		z[i] = 0.0;
		I[i] = 0.0;
	}

	/* for a feasible basis, compute the solution */
	if ( Result == LCP_FEASIBLE || Result == LCP_PRETERMINATED )
	{
#ifdef MATLAB_MEX_FILE
		t1 = clock();
#endif
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
			I[j] = Index_Get(pB->pZ, i)+m+1;
		}
#ifdef MATLAB_MEX_FILE
		t2 = clock();
		total_time+=(double)(t2-t1)/(double)CLOCKS_PER_SEC;
		if (options.clock) {
			printf("Time in total needed to solve LCP: %lf seconds\n",
			       total_time + (double)(t2-t1)/(double)CLOCKS_PER_SEC);
		}
#endif
		
		if (options.normalize) {
			/* do the backward normalization */
			/* z = diag(c)*zn */
			for (i=0; i<m; i++) {
				z[i] = c[i]*z[i];
			}
			
			/* since the normalization does not compute w properly, we recalculate it from
			 * recovered z */
			/* write data to Mn = M */
			memcpy(Mn, *M, n*sizeof(double));
			/* write data to qn = q */
			memcpy(qn, *q, m*sizeof(double));
			/* copy w <- q; */
			dcopy(&m, qn, &inc, w, &inc);
			/* compute w = M*z + q */
			dgemv(&T, &m, &m, &alpha, Mn, &m, z, &inc, &alpha, w, &inc);
			/* if w is less than eps, consider it as zero */
			for (i=0; i<m; i++) {
				if (w[i]<options.zerotol) {
					w[i] = 0.0;
				}
			}
		}	       
	}


	/* outputs */
	*exitflag = (real_T )Result;
	*pivots =(real_T)pivs;
	*time = (real_T)total_time;

	/* reset dimensions and values in basis for a recursive call */
	Reinitialize_Basis(m, pB);
	
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


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
