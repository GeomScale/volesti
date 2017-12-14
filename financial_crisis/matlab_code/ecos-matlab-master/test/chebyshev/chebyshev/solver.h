/* stuff about open source license
 * ....
 * The problem specification for this solver is: 
 *
 *     
 *     parameter a1t(1,2)
 *     parameter a2t(1,2)
 *     parameter a3t(1,2)
 *     parameter a4t(1,2)
 *     
 *     parameter na1
 *     parameter na2
 *     parameter na3
 *     parameter na4
 *     
 *     parameter b1
 *     parameter b2
 *     parameter b3
 *     parameter b4
 *     
 *     variable r
 *     variable x(2)
 *     
 *     maximize r
 *     subject to
 *         a1t*x + na1*r <= b1
 *         a2t*x + na2*r <= b2
 *         a3t*x + na3*r <= b3
 *         a4t*x + na4*r <= b4
 *
 * For now, parameters are *dense*, and we don't respect sparsity. The sparsity
 * structure has to be specified during code generation time. We could do the more
 * generic thing and allow sparse matrices, but as a first cut, we won't.
 * Version 0.0.1
 * Eric Chu, Alex Domahidi, Neal Parikh, Stephen Boyd (c) 2012 or something...
 */

#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "ecos.h"

/*
 * struct params_t (or `params`)
 * =============================
 * This structure contains the data for all parameters in your problem.
 *
 */
typedef struct params_t {
  double a1t[1][2];
  double a2t[1][2];
  double a3t[1][2];
  double a4t[1][2];
  double b1;
  double b2;
  double b3;
  double b4;
  double na1;
  double na2;
  double na3;
  double na4;
} params;

/*
 * struct vars_t (or `vars`)
 * =========================
 * This structure stores the solution variables for your problem.
 *
 */
typedef struct vars_t {
  double r;
  double x[2];
} vars;

pwork *setup(params *p);        /* setting up workspace (assumes params already declared) */
int solve(pwork *w, vars *sol); /* solve the problem (assumes vars already declared) */
void cleanup(pwork *w);         /* clean up workspace */

#endif    /* solver.h */
