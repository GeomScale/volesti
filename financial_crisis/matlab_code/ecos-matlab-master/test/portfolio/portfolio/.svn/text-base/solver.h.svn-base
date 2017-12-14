/* stuff about open source license
 * ....
 * The problem specification for this solver is: 
 *
 *     dimension n = 10  # dimension of problem
 *     dimension m = 3
 *     
 *     # testing comments
 *     variable x(n)
 *     
 *     parameter mu(n)
 *     parameter B
 *     parameter gamma positive
 *     parameter F(n,m)
 *     parameter d(n)
 *     
 *     # ignoring comments?
 *     maximize mu'*x - gamma*(square(norm(F'*x)) + square(norm(diag(d)*x)))
 *     subject to
 *         sum(x) == B # ignored?
 *         x >= 0
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
  double B;
  double F[10][3];
  double d[10];
  double gamma;
  double mu[10];
} params;

/*
 * struct vars_t (or `vars`)
 * =========================
 * This structure stores the solution variables for your problem.
 *
 */
typedef struct vars_t {
  double x[10];
} vars;

pwork *setup(params *p);        /* setting up workspace (assumes params already declared) */
int solve(pwork *w, vars *sol); /* solve the problem (assumes vars already declared) */
void cleanup(pwork *w);         /* clean up workspace */

#endif    /* solver.h */
