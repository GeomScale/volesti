/* stuff about open source license
 * ....
 * The problem specification for this solver is: 
 *
 *     dimension m = 10
 *     dimension n = 2
 *     
 *     parameter gamma positive
 *     
 *     parameter X(m,n)
 *     parameter Y(m,n)
 *     
 *     variable a(2)
 *     variable b
 *     
 *     minimize (norm(a) + gamma*sum(pos(1 - X*a + b) + pos(1 + Y*a - b)))
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
  double X[10][2];
  double Y[10][2];
  double gamma;
} params;

/*
 * struct vars_t (or `vars`)
 * =========================
 * This structure stores the solution variables for your problem.
 *
 */
typedef struct vars_t {
  double a[2];
  double b;
} vars;

pwork *setup(params *p);        /* setting up workspace (assumes params already declared) */
int solve(pwork *w, vars *sol); /* solve the problem (assumes vars already declared) */
void cleanup(pwork *w);         /* clean up workspace */

#endif    /* solver.h */
