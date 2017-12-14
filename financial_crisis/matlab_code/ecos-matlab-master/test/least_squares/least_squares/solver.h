/* stuff about open source license
 * ....
 * The problem specification for this solver is: 
 *
 *     dimension m = 10
 *     dimension n = 5
 *     dimension p = 3
 *     
 *     parameter A(m,n)
 *     parameter b(m)
 *     parameter C(p,n)
 *     parameter d(p)
 *     parameter e
 *     
 *     variable x(n)
 *     
 *     minimize square(norm(A*x - b))
 *     subject to
 *         C*x == d
 *         norm_inf(x) <= e
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
  double A[10][5];
  double C[3][5];
  double b[10];
  double d[3];
  double e;
} params;

/*
 * struct vars_t (or `vars`)
 * =========================
 * This structure stores the solution variables for your problem.
 *
 */
typedef struct vars_t {
  double x[5];
} vars;

pwork *setup(params *p);        /* setting up workspace (assumes params already declared) */
int solve(pwork *w, vars *sol); /* solve the problem (assumes vars already declared) */
void cleanup(pwork *w);         /* clean up workspace */

#endif    /* solver.h */
