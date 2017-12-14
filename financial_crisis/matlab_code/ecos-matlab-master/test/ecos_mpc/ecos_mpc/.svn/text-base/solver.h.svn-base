/* stuff about open source license
 * ....
 * The problem specification for this solver is: 
 *
 *     variable x1
 *     variable x2
 *     variable x3
 *     variable x4
 *     variable x5
 *     variable x6
 *     variable x7
 *     
 *     parameter p positive
 *     parameter q positive
 *     parameter r positive
 *     parameter a
 *     parameter b
 *     parameter x_0
 *     parameter one_over_q
 *     
 *     minimize (q*square(x1)  + r*square(x5) + q*max([x2;1]) - r*min([x6;0]) + quad_over_lin(x3, one_over_q) - r*sqrt(x7) + p*inv_pos(x4))
 *     subject to
 *     x2 == a*x1 + b*x5
 *     x3 == a*x2 + b*x6
 *     x4 == a*x3 + b*x7
 *     x1 == x_0
 *     0 <= x5
 *     x5 <= 1
 *     0 <= x6
 *     x6 <= 1
 *     0 <= x7
 *     x7 <= 1
 *     
 *     
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
  double a;
  double b;
  double one_over_q;
  double p;
  double q;
  double r;
  double x_0;
} params;

/*
 * struct vars_t (or `vars`)
 * =========================
 * This structure stores the solution variables for your problem.
 *
 */
typedef struct vars_t {
  double x1;
  double x2;
  double x3;
  double x4;
  double x5;
  double x6;
  double x7;
} vars;

pwork *setup(params *p);        /* setting up workspace (assumes params already declared) */
int solve(pwork *w, vars *sol); /* solve the problem (assumes vars already declared) */
void cleanup(pwork *w);         /* clean up workspace */

#endif    /* solver.h */
