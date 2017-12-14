/* stuff about open source license
 * ....
 * The problem specification for this solver is: 
 *
 *     variable x1
 *     variable x2
 *     
 *     parameter a11
 *     parameter a12
 *     parameter a21
 *     parameter a22
 *     parameter a31
 *     parameter a32
 *     parameter a41
 *     parameter a42
 *     
 *     parameter b1
 *     parameter b2
 *     parameter b3
 *     parameter b4
 *     
 *     parameter c11
 *     parameter c12
 *     parameter d
 *     parameter e
 *     
 *     minimize square(a11*x1 + a12*x2 - b1) + square(a21*x1 + a22*x2 - b2) + square(a31*x1 + a32*x2 - b3) + square(a41*x1 + a42*x2 - b4)
 *     subject to
 *         c11*x1 + c12*x2 == d
 *         max([abs(x1);abs(x2)]) <= e
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
  double a11;
  double a12;
  double a21;
  double a22;
  double a31;
  double a32;
  double a41;
  double a42;
  double b1;
  double b2;
  double b3;
  double b4;
  double c11;
  double c12;
  double d;
  double e;
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
} vars;

pwork *setup(params *p);        /* setting up workspace (assumes params already declared) */
int solve(pwork *w, vars *sol); /* solve the problem (assumes vars already declared) */
void cleanup(pwork *w);         /* clean up workspace */

#endif    /* solver.h */
