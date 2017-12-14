/* cddmex.c: Main program to call cdd library from Matlab (as a mex file).
   For details of proper Matlab call check accompanying cddmex.m file.

  revised 2010 by M. Herceg, ETH Zurich
  (C) 2004 by M. Kvasnica, ETH Zurich
  (C) 2003 by M. Baotic, Zurich, October 14, 2003
  (C) 2002 by F. Torrisi and M. Baotic, Zurich, 2002
*/

/* The cdd library cddlib-094g was written by
   Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.94g, 03/30/2012
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
*/

/*  This program is free software; you can redistribute it and/or modify
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


#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mex.h"

#define CDDMEX_VERSION "1.0"

mxArray * MB_set_LPsol_MatrixPtr(const dd_LPPtr lp)
{
	mxArray * P;	
	mxArray * tmpx;
	mxArray * tmpy;	
	mxArray * tmpobj;
	mxArray * tmphow;
	mxArray * tmpactive;
	double * x;
	double * y;
	double * obj;
	double * how;
	double * active;
	int ii, jj;
	int nactive;
	

 	if (lp !=NULL) {
    		const char *f[] = {"xopt", "lambda", "how", "objlp", "active"};
    		mwSize dims[] = {1};

    		P = mxCreateStructArray(1, dims, 5, f);
    		
    		/* store primal solution */
    		tmpx = mxCreateDoubleMatrix((lp->d)-1, 1, mxREAL);
    		x = mxGetPr(tmpx);
    		for (ii=1; ii<=(lp->d)-1; ii++)
    			x[ii-1]=(double)(lp->sol[ii])[0];
		mxSetField(P, 0, "xopt", tmpx);

    		/* store Lagrange multipliers (i.e. dual solution) */
    		tmpy = mxCreateDoubleMatrix((lp->m)-1, 1, mxREAL);
    		y = mxGetPr(tmpy);
    		for (ii=1; ii<lp->m; ii++)
    			y[ii-1]=0;
    		nactive=0;
		for (ii=1; ii<lp->d; ii++){
			if (lp->nbindex[ii+1]>0) {
				nactive++;
				y[lp->nbindex[ii+1]-1]=(lp->dsol[ii])[0];
			}
		}
		mxSetField(P, 0, "lambda", tmpy);

		/* store objective */
		tmphow = mxCreateDoubleMatrix(1, 1, mxREAL);
    		how = mxGetPr(tmphow);
		how[0]=lp->LPS;
      		mxSetField(P, 0, "how", tmphow);

		/* store objective value */
		tmpobj = mxCreateDoubleMatrix(1, 1, mxREAL);
    		obj = mxGetPr(tmpobj);
		obj[0]=(lp->optvalue)[0];
      		mxSetField(P, 0, "objlp", tmpobj);

		/* store nonbasis */
    		tmpactive = mxCreateDoubleMatrix(nactive, 1, mxREAL);
    		active = mxGetPr(tmpactive);
    		jj=0;
		for (ii=1; ii<lp->d; ii++){
			if (lp->nbindex[ii+1]>0){
				active[jj]=(double)lp->nbindex[ii+1];
				jj++;
			}
		}
		mxSetField(P, 0, "active", tmpactive);


      		return P;
	}
	return 0;
}

	
dd_MatrixPtr FT_get_H_MatrixPtr(const mxArray * in)
{
	mxArray * tmpa;	
	mxArray * tmpb;	
	mxArray * tmpl;	
	dd_MatrixPtr A;
	int eqsize, m, n, i, j;
	double * a;
	double * b;
	double * lin;
	
	if ((tmpa = mxGetField(in, 0, "A")) && 
	    (tmpb = mxGetField(in, 0, "B")) &&
	    (mxGetNumberOfDimensions(tmpa) <= 2) &&
	    (mxGetNumberOfDimensions(tmpb) <= 2) &&
	    (mxGetM(tmpa) == mxGetM(tmpb)) &&
	    (mxGetN(tmpb) == 1)) {
		m = mxGetM(tmpa);
		n = mxGetN(tmpa) + 1;
		a = mxGetPr(tmpa);
		b = mxGetPr(tmpb);		
		A = dd_CreateMatrix(m, n);
		for (i = 0; i < m; i++) {
			dd_set_d(A->matrix[i][0],b[i]);
			for (j = 0; j < n - 1; j++) {
				dd_set_d(A->matrix[i][j + 1],- a[j * m + i]);
			}
		}
		A->representation = dd_Inequality;		
		A->numbtype = dd_Real;
		/* set lineality if present */
		if ((tmpl = mxGetField(in, 0, "lin")) &&
		    (mxGetNumberOfDimensions(tmpl) <= 2) && 
		    (mxGetM(tmpl) == 1)) {
		    	lin = mxGetPr(tmpl);
		    	eqsize = mxGetN(tmpl);
		    	for (i = 0; i < eqsize; i++) {
		    		if (lin[i] <= (double) m) 
     					set_addelem(A->linset,(long) lin[i]);
     				else
     					mexWarnMsgTxt("Error in the lineality vector ");
  			}
  		}
		return A;
	}
	return 0;
}


dd_LPPtr MB_get_LP_MatrixPtr(const mxArray * in)
{
	mxArray * tmpobj;
	dd_MatrixPtr A;
	int j;
	double * value;
	dd_ErrorType error=dd_NoError;
	dd_LPPtr lp;
	
	dd_set_global_constants(); /* First, this must be called once to use cddlib. */
	
	if (A=FT_get_H_MatrixPtr(in)) {
		/* set objective */
		if ((tmpobj = mxGetField(in, 0, "obj")) &&
		    (mxGetNumberOfDimensions(tmpobj) <= 2) && 
		    (mxGetM(tmpobj) == 1)) {
		    	value = mxGetPr(tmpobj);
  			A->objective = dd_LPmin;
  			dd_set_d(A->rowvec[0], 0.0); /* there is no constant term */
  			for (j = 1; j < A->colsize; j++)
  				dd_set_d(A->rowvec[j],value[j-1]);
  			lp=dd_Matrix2LP(A, &error);
  			dd_FreeMatrix(A);
			return lp;
  		}else{
  			mexErrMsgTxt("Error in the setting of LP objective.");
  		}
	}else{
		mexErrMsgTxt("Error in the setting of LP matrix.");
	}
	return 0;
}

mxArray * FT_set_V_MatrixPtr(const dd_MatrixPtr M)
{
	mxArray * P;
	int mr, mv, m, i, j;
	int ir, iv;
	mxArray * tmpv;	
	mxArray * tmpr;	
	mxArray * tmpvpos;	
	mxArray * tmprpos;	
	mxArray * tmpl;	
	
	double * v;
	double * r; 
	double * vpos; 
	double * rpos; 
	double * l; 

	int k=0;
	dd_rowrange ii;
		
 	if ((M !=NULL) &&
 	    (M->representation == dd_Generator)) {
    		const char *f[] = {"V", "R", "rpos", "vpos", "lin"};
    		mwSize dims[] = {1};

    		P = mxCreateStructArray(1, dims, 5, f);
		if (set_card(M->linset)) {
			tmpl = mxCreateDoubleMatrix(1, set_card(M->linset), mxREAL);
			l = mxGetPr(tmpl);		
    			for (ii=1; ii<=M->rowsize; ii++) {
      				if (set_member(ii, M->linset)) {
      					l[k] = (double)ii;
      					k ++;
      				}
      			}
      			mxSetField(P, 0, "lin", tmpl);
		}
		/* Count the rays and vertices */
		mr = 0;
		mv = 0;
		for (i = 0 ; i < (int)(M->rowsize); i++) {
			if (dd_get_d(M->matrix[i][0]) == 0) {
				mr++;
			} else { 
				mv++;
			}
		}
		m = mr + mv;
		/* Allocate the space in MATLAB */
		tmpr = mxCreateDoubleMatrix(mr, M->colsize - 1, mxREAL);
		r = mxGetPr(tmpr);
		tmpv = mxCreateDoubleMatrix(mv, M->colsize - 1, mxREAL);
		v = mxGetPr(tmpv);
		tmprpos = mxCreateDoubleMatrix(mr, 1, mxREAL);
		rpos = mxGetPr(tmprpos);
		tmpvpos = mxCreateDoubleMatrix(mv, 1, mxREAL);
		vpos = mxGetPr(tmpvpos);
		ir = 0; iv = 0;
	  	for (i = 0 ; i < m; i++) {
			if (dd_get_d(M->matrix[i][0]) == 0) {
				/* This is a ray */
    				for (j = 0; j < (int)(M->colsize) - 1; j++) {
      					r[ir + j * mr] = dd_get_d(M->matrix[i][j + 1]);
      				}
      				rpos[ir] = i + 1;
				ir++;
			} else { 
				/* This is a vertex*/
    				for (j = 0; j < (int)(M->colsize) - 1; j++) {
      					v[iv + j * mv] = dd_get_d(M->matrix[i][j + 1]);
      				}
      				vpos[iv] = i + 1;
				iv++;
			}

      		}		
      		mxSetField(P, 0, "V", tmpv);
      		mxSetField(P, 0, "R", tmpr);
      		mxSetField(P, 0, "vpos", tmpvpos);
      		mxSetField(P, 0, "rpos", tmprpos);
      		return P;
	}
	return 0;
}


dd_MatrixPtr FT_get_V_MatrixPtr(const mxArray * in)
{
	mxArray * tmpv;	
	mxArray * tmpr;	
	dd_MatrixPtr V;
	int mr, m, n, i, j;
	double * v;
	double * r;
	
	if ((tmpv = mxGetField(in, 0, "V")) && 
	    (mxGetNumberOfDimensions(tmpv) <= 2)) {
	    	
		if ((tmpr = mxGetField(in, 0, "R")) && 
	    	    (mxGetNumberOfDimensions(tmpv) <= 2)) {
	    	    	mr = mxGetM(tmpr);
	    	    	r = mxGetPr(tmpr);
		} else mr = 0;
	    		    
		m = mxGetM(tmpv);
		n = mxGetN(tmpv) + 1;
		v = mxGetPr(tmpv);
		V = dd_CreateMatrix(m + mr, n);
		for (i = 0; i < m; i++) {
			dd_set_si(V->matrix[i][0],1);
			for (j = 0; j < n - 1; j++) {
				dd_set_d(V->matrix[i][j + 1], v[i + j * m]);
			}
		}
		for (i = m; i < m + mr; i++) {
			dd_set_si(V->matrix[i][0],0);
			for (j = 0; j < n - 1; j++) {
				dd_set_d(V->matrix[i][j + 1], r[(i - m) + j * mr]);
			}
		}		
		V->representation = dd_Generator;	
		V->numbtype = dd_Real;
		return V;
	}
	return 0;
}

mxArray * FT_set_SetFamilyPtr(const dd_SetFamilyPtr F)
{
	mxArray * A, * tmp;
	double * r;
	int i, j, k;
	int num;
	
	if (F) {
		num = (int)(F->famsize);
		A = mxCreateCellMatrix(num, 1); 
		for (i=0; i<(int)(F->famsize); i++) {
			num = (int)set_card(F->set[i]);
			tmp = mxCreateDoubleMatrix(1, num, mxREAL);
			r = mxGetPr(tmp);
			k = 0;
			for (j = 1; j <= *(F->set[i]); j++) {
				if (set_member(j,F->set[i]))
					r[k++] = j;
			}
			mxSetCell(A, i, tmp);
		}
		return A;
	}
	return 0;
}

mxArray * FT_set_Set(const set_type S)
{
	mxArray * tmp;
	double * r;
	int i, j;
	int num;
	
	if (S) {
		num = (int)(S[0]);
		j = 0;
		for (i = 1; i <= num; i++) {
			if (set_member(i,S)) {
				j++; 
			}
		}
		tmp = mxCreateDoubleMatrix(1, j, mxREAL);
		r = mxGetPr(tmp);
		j = 0;
		for (i = 1; i <= num; i++) {
			if (set_member(i,S)) {
				r[j++] = i;
			}
		}	
		return tmp;
	}
	return 0;
}

mxArray * FT_set_H_MatrixPtr(const dd_MatrixPtr M)
{
	mxArray * P;	
	mxArray * tmpa;	
	mxArray * tmpb;	
	mxArray * tmpl;
	int i, j;
	double * a;
	double * b;
	double * l;
	int k=0;

	dd_rowrange ii;

 	if ((M !=NULL) &&
 	    (M->representation == dd_Inequality)) {
    		const char *f[] = {"A", "B", "lin"};
    		mwSize dims[] = {1};

    		P = mxCreateStructArray(1, dims, 3, f);
		if (set_card(M->linset)) {
			tmpl = mxCreateDoubleMatrix(1, set_card(M->linset), mxREAL);
			l = mxGetPr(tmpl);		
    			for (ii=1; ii<=M->rowsize; ii++) {
      				if (set_member(ii, M->linset)) {
      					l[k] = (int)ii;
      					k ++;
      				}
      			}
      			mxSetField(P, 0, "lin", tmpl);
		}

		tmpb = mxCreateDoubleMatrix(M->rowsize, 1, mxREAL);
		b = mxGetPr(tmpb);
		tmpa = mxCreateDoubleMatrix(M->rowsize, M->colsize - 1, mxREAL);
		a = mxGetPr(tmpa);
	  	for (i = 0 ; i < (int)(M->rowsize); i++) {
	  		b[i] = dd_get_d(M->matrix[i][0]);
    			for (j = 0; j < (int)(M->colsize) - 1; j++) {
      				a[i + j * (int)(M->rowsize)] = - dd_get_d(M->matrix[i][j + 1]);
      			}
      		}
      		mxSetField(P, 0, "A", tmpa);
      		mxSetField(P, 0, "B", tmpb);
      		return P;
	}
	return 0;
}

void
hull(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_PolyhedraPtr P;
	dd_ErrorType err;
	dd_MatrixPtr H,V;
	
	if (nrhs  == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
		V = FT_get_V_MatrixPtr(prhs[0]);		
		dd_set_global_constants();  /* First, this must be called. */

		P = dd_DDMatrix2Poly(V, &err); /* compute the second representation */
		if (err == dd_NoError) {
			H = dd_CopyInequalities(P);
			plhs[0] = FT_set_H_MatrixPtr(H);
			dd_FreeMatrix(H);
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");
  		}
		dd_FreeMatrix(V);
  		dd_FreePolyhedra(P);
  		return;
	} else {
		mexErrMsgTxt("hull expects a V input struct and produces an H output struct");
	}
}

void
v_hull_extreme(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_PolyhedraPtr P,P1;
	dd_ErrorType err;
	dd_MatrixPtr H,V,V1;
	
	if (nrhs  == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
		V = FT_get_V_MatrixPtr(prhs[0]);		
		dd_set_global_constants();  /* First, this must be called. */

		P = dd_DDMatrix2Poly(V, &err); /* compute the second representation */
		if (err == dd_NoError) {
			H = dd_CopyInequalities(P);
			P1 = dd_DDMatrix2Poly(H, &err); /* compute the second representation */
			if (err == dd_NoError) {
				V1 = dd_CopyGenerators(P1);
				plhs[0] = FT_set_V_MatrixPtr(V1);
				dd_FreeMatrix(V1);
			} else {
    				dd_WriteErrorMessages(stdout,err);
    				mexErrMsgTxt("CDD returned an error, see above(!) for details");    			
  			}
			dd_FreePolyhedra(P1);
			dd_FreeMatrix(H);
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");
  		}
		dd_FreeMatrix(V);
  		dd_FreePolyhedra(P);
  		return;
	} else {
		mexErrMsgTxt("v_hull_extreme expects a V input struct and produces a V output struct");
	}
}

void
extreme(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_PolyhedraPtr P;
	dd_ErrorType err;
	dd_MatrixPtr H,V;
	
	if (nrhs  == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
		dd_set_global_constants();  /* First, this must be called. */
		H = FT_get_H_MatrixPtr(prhs[0]);		

		P = dd_DDMatrix2Poly(H, &err); /* compute the second representation */
		if (err == dd_NoError) {
			V = dd_CopyGenerators(P);
			plhs[0] = FT_set_V_MatrixPtr(V);
			dd_FreeMatrix(V);
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");    			
  		}
		dd_FreeMatrix(H);
  		dd_FreePolyhedra(P);
  		return;
	} else {
		mexErrMsgTxt("extreme expects an H input struct and produces a V output struct");
	}
}

void
adj_extreme(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_PolyhedraPtr P;
	dd_ErrorType err;
	dd_MatrixPtr H,V;
	dd_SetFamilyPtr A;
	
	if (nrhs  == 1 && nlhs == 2 && mxIsStruct(prhs[0])) {
		dd_set_global_constants();  /* First, this must be called. */
		H = FT_get_H_MatrixPtr(prhs[0]);		

		P = dd_DDMatrix2Poly(H, &err); /* compute the second representation */
		if (err == dd_NoError) {
			V = dd_CopyGenerators(P);
			A = dd_CopyAdjacency(P);
			plhs[0] = FT_set_V_MatrixPtr(V);
			plhs[1] = FT_set_SetFamilyPtr(A);
			dd_FreeMatrix(V);
			dd_FreeSetFamily(A);
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");    			
  		}
		dd_FreeMatrix(H);
  		dd_FreePolyhedra(P);
  		return;
	} else {
		mexErrMsgTxt("adj_extreme expects an H input struct and produces a V output struct and the adjacency struct");
	}
}

void
reduce_h(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_ErrorType err;
	dd_MatrixPtr H,H1;
	dd_rowset red;
	
	if (nrhs  == 1 && nlhs >= 1 && nlhs <= 2 && mxIsStruct(prhs[0])) {
		dd_set_global_constants();  /* First, this must be called. */
		H = FT_get_H_MatrixPtr(prhs[0]);		
		red = dd_RedundantRows(H, &err); /* find redundant rows */
		if (err == dd_NoError) {
			/* remove the red rows */
			H1 = dd_MatrixSubmatrix(H, red);
			plhs[0] = FT_set_H_MatrixPtr(H1);
			dd_FreeMatrix(H1);
			if (nlhs == 2) {
				plhs[1] = FT_set_Set(red);
			}
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");
  		}
		dd_FreeMatrix(H);
		set_free(red);
  		return;
	} else {
		mexErrMsgTxt("reduce_h expects an H input struct and produces a H output struct and an optional vector of removed rows");
	}
}

void
reduce_v(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_ErrorType err;
	dd_MatrixPtr V,V1;
	dd_rowset red;
	
	if (nrhs  == 1 && nlhs >= 1 && nlhs <= 2 && mxIsStruct(prhs[0])) {
		dd_set_global_constants();  /* First, this must be called. */
		V = FT_get_V_MatrixPtr(prhs[0]);		
		red = dd_RedundantRows(V, &err); /* find redundant rows */
		if (err == dd_NoError) {
			/* remove the red rows */
			V1 = dd_MatrixSubmatrix(V, red);
			plhs[0] = FT_set_V_MatrixPtr(V1);
			dd_FreeMatrix(V1);
			if (nlhs == 2) {
				plhs[1] = FT_set_Set(red);
			}
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");
  		}
		dd_FreeMatrix(V);
		set_free(red);
  		return;
	} else {
		mexErrMsgTxt("reduce_v expects an V input struct and produces a V output struct and an optional vector of removed vertices");
	}
}

void
copy_v(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_MatrixPtr V;
	
	if (nrhs  == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
		dd_set_global_constants();  /* First, this must be called. */
		V = FT_get_V_MatrixPtr(prhs[0]);		
		plhs[0] = FT_set_V_MatrixPtr(V);
		dd_FreeMatrix(V);
  		return;
	} else {
		mexErrMsgTxt("copy_v expects a V input struct and produces a V output struct");
	}
}

void
copy_h(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_MatrixPtr H;
	
	if (nrhs  == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
		dd_set_global_constants();  /* First, this must be called. */
		H = FT_get_H_MatrixPtr(prhs[0]);		
		plhs[0] = FT_set_H_MatrixPtr(H);
		dd_FreeMatrix(H);
  		return;
	} else {
		mexErrMsgTxt("copy_h expects a H input struct and produces a H output struct");
	}
}

void
solve_lp(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	/* The original LP data  m x n matrix 
	   = | b   -A  |
	     | c0  c^T |,

	where the LP to be solved is to
	maximize  c^T x  +   c0
	subj. to
	     A   x  <=  b.
	*/

	dd_ErrorType error=dd_NoError;
	dd_LPSolverType solver=dd_CrissCross; /* either DualSimplex or CrissCross */
	dd_LPPtr lp;   /* pointer to LP data structure that is not visible by user. */
  
	dd_set_global_constants(); /* First, this must be called once to use cddlib. */


	/* Input an LP using the cdd library  */
	lp = MB_get_LP_MatrixPtr(prhs[0]);

    
	/* Solve the LP by cdd LP solver. */
	dd_LPSolve(lp,solver,&error);
	if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
	
	/* Take the solution. */
	plhs[0] = MB_set_LPsol_MatrixPtr(lp);


	/* Free allocated spaces. */
	dd_FreeLPData(lp);
}

void
solve_lp_DS(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  /* uses Dual Simplex method */
	/* The original LP data  m x n matrix 
	   = | b   -A  |
	     | c0  c^T |,

	where the LP to be solved is to
	maximize  c^T x  +   c0
	subj. to
	     A   x  <=  b.
	*/
  
	dd_ErrorType error=dd_NoError;
	dd_LPSolverType solver=dd_DualSimplex;
	dd_LPPtr lp;   /* pointer to LP data structure that is not visible by user. */
  
	dd_set_global_constants(); /* First, this must be called once to use cddlib. */


	/* Input an LP using the cdd library  */
	lp = MB_get_LP_MatrixPtr(prhs[0]);

    
	/* Solve the LP by cdd LP solver. */
	dd_LPSolve(lp,solver,&error);
	if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
	
	/* Take the solution. */
	plhs[0] = MB_set_LPsol_MatrixPtr(lp);


	/* Free allocated spaces. */
	dd_FreeLPData(lp);
}


void
find_interior(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	/* We would like to find an iterior point
	   for a polyhedron in H representation 
	     A   x  <=  b.
	*/

	dd_ErrorType error=dd_NoError;
	dd_LPSolverType solver=dd_CrissCross; /* either DualSimplex or CrissCross */
	dd_LPPtr lp, lp1;   /* pointer to LP data structure that is not visible by user. */
	dd_MatrixPtr A;
	int j;  
  
	dd_set_global_constants(); /* First, this must be called once to use cddlib. */


	/* Input an LP using the cdd library  */
	/* lp = MB_get_LP_MatrixPtr(prhs[0]); */

	if (A=FT_get_H_MatrixPtr(prhs[0])) {
		/* set objective */
		A->objective = dd_LPmin;
		for (j = 0; j < A->colsize; j++)
			dd_set_d(A->rowvec[j],0.0);
		lp=dd_Matrix2LP(A, &error);
  		dd_FreeMatrix(A);
	}else{
		mexErrMsgTxt("Error in the setting of LP matrix.");
	}

	
	lp1=dd_MakeLPforInteriorFinding(lp);
	dd_LPSolve(lp1,solver,&error);
	if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
	    
	/* Take the solution. */
	plhs[0] = MB_set_LPsol_MatrixPtr(lp1);


	/* Free allocated spaces. */
	dd_FreeLPData(lp);
	dd_FreeLPData(lp1);
}


void
find_interior_DS(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  /* uses Dual Simplex method */
	/* We would like to find an iterior point
	   for a polyhedron in H representation 
	     A   x  <=  b.
	*/

	dd_ErrorType error=dd_NoError;
	dd_LPSolverType solver=dd_DualSimplex;
	dd_LPPtr lp, lp1;   /* pointer to LP data structure that is not visible by user. */
	dd_MatrixPtr A;
	int j;  
  
	dd_set_global_constants(); /* First, this must be called once to use cddlib. */


	/* Input an LP using the cdd library  */
	/* lp = MB_get_LP_MatrixPtr(prhs[0]); */

	if (A=FT_get_H_MatrixPtr(prhs[0])) {
		/* set objective */
		A->objective = dd_LPmin;
		for (j = 0; j < A->colsize; j++)
			dd_set_d(A->rowvec[j],0.0);
		lp=dd_Matrix2LP(A, &error);
  		dd_FreeMatrix(A);
	}else{
		mexErrMsgTxt("Error in the setting of LP matrix.");
	}

	
	lp1=dd_MakeLPforInteriorFinding(lp);
	dd_LPSolve(lp1,solver,&error);
	if (error!=dd_NoError) dd_WriteErrorMessages(stdout, error);
	    
	/* Take the solution. */
	plhs[0] = MB_set_LPsol_MatrixPtr(lp1);


	/* Free allocated spaces. */
	dd_FreeLPData(lp);
	dd_FreeLPData(lp1);
}

void
mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int buflen, status;
	char *input_buf;
	
	if (nrhs<1){
		mexErrMsgTxt("Input 1 must be a row vector string");
	}
		
	
	/* 1. input must be a string and row vector */
	if (mxIsChar(prhs[0]) && (mxGetM(prhs[0]) == 1)) {
	/* get the length of the input string */
		buflen = mxGetN(prhs[0]) + 1;
		/* allocate memory for input and output strings */
    		input_buf=mxCalloc(buflen, sizeof(char));
		/* copy the string data from prhs[0] into a C string input_ buf. */
    		status = mxGetString(prhs[0], input_buf, buflen);
		if (strcmp(input_buf,"hull\0") == 0) {
			hull(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"extreme\0") == 0) {
			extreme(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"reduce_h\0") == 0) {
			reduce_h(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"reduce_v\0") == 0) {
			reduce_v(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"copy_v\0") == 0) {
			copy_v(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"copy_h\0") == 0) {
			copy_h(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"v_hull_extreme\0") == 0) {
			v_hull_extreme(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"adj_extreme\0") == 0) {
			adj_extreme(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"solve_lp\0") == 0) {
			solve_lp(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"solve_lp_DS\0") == 0) {
			solve_lp_DS(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"find_interior\0") == 0) {
			find_interior(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"find_interior_DS\0") == 0) {
			find_interior_DS(nlhs, plhs, nrhs -1, prhs + 1);
			return;
		}
		if (strcmp(input_buf,"version\0") == 0) {
			printf("Version %s\n", CDDMEX_VERSION);
			return;
		}
		mexErrMsgTxt("Unknown function");		
	} else {
		mexErrMsgTxt("Input 1 must be a row vector string");
	}	
}

void
adjacency(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	dd_PolyhedraPtr P;
	dd_ErrorType err;
	dd_MatrixPtr V;
	dd_SetFamilyPtr A;
	
		
	if (mxIsStruct(prhs[0])) {
		V = FT_get_V_MatrixPtr(prhs[0]);		
		dd_set_global_constants();  /* First, this must be called. */

		P = dd_DDMatrix2Poly(V, &err); /* compute the second representation */
		if (err == dd_NoError) {
			A = dd_CopyInputAdjacency(P);
			plhs[0] = FT_set_SetFamilyPtr(A);
			dd_FreeSetFamily(A);
		} else {
    			dd_WriteErrorMessages(stdout,err);
    			mexErrMsgTxt("CDD returned an error, see above(!) for details");
  		}
		dd_FreeMatrix(V);
  		dd_FreePolyhedra(P);
	}
}
