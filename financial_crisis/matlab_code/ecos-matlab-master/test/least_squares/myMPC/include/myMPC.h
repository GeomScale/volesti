/*
FORCES - Fast interior point code generation for multistage problems.
Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
Automatic Control Laboratory, ETH Zurich.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __myMPC_H__
#define __myMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double myMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef myMPC_SET_PRINTLEVEL
#define myMPC_SET_PRINTLEVEL    (2)
#endif

/* maximum number of iterations  */
#define myMPC_SET_MAXIT         (30)	

/* scaling factor of line search (affine direction) */
#define myMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define myMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define myMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define myMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define myMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define myMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define myMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define myMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define myMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define myMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define myMPC_NOPROGRESS   (-7)


/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct myMPC_params
{
    /* vector of size 4 */
    myMPC_FLOAT z1[4];

} myMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct myMPC_output
{
    /* vector of size 1 */
    myMPC_FLOAT u1[1];

} myMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct myMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    myMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    myMPC_FLOAT res_ineq;

    /* primal objective */
    myMPC_FLOAT pobj;	
	
    /* dual objective */
    myMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    myMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    myMPC_FLOAT rdgap;		

    /* duality measure */
    myMPC_FLOAT mu;

	/* duality measure (after affine step) */
    myMPC_FLOAT mu_aff;
	
    /* centering parameter */
    myMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    myMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    myMPC_FLOAT step_cc;    

	/* solvertime */
	myMPC_FLOAT solvetime;   

} myMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int myMPC_solve(myMPC_params* params, myMPC_output* output, myMPC_info* info);


#endif