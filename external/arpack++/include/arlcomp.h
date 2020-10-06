/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLComp.h.
   ALTERED copy of dcomplex.h and scomplex.h (from SuperLU package).
   Structure complex was renamed to lscomplex.
   Structure doublecomplex was renamed to ldcomplex.
*/

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 * </pre>
 */
#ifndef __SUPERLU_DCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_DCOMPLEX

/* 
 * This header file is to be included in source files z*.c
 */
#ifndef DCOMPLEX_INCLUDE
#define DCOMPLEX_INCLUDE

typedef struct { double r, i; } ldcomplex;


/* Macro definitions */

/*! \brief Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/*! \brief Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/*! \brief Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/*! \brief Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

#define zz_conj(a, b) { \
        (a)->r = (b)->r; \
        (a)->i = -((b)->i); \
    }

/*! \brief Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in dcomplex.c */
void z_div(ldcomplex *, ldcomplex *, ldcomplex *);
double z_abs(ldcomplex *);     /* exact */
double z_abs1(ldcomplex *);    /* approximate */
void z_exp(ldcomplex *, ldcomplex *);
void d_cnjg(ldcomplex *r, ldcomplex *z);
double d_imag(ldcomplex *);
ldcomplex z_sgn(ldcomplex *);
ldcomplex z_sqrt(ldcomplex *);


#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_DCOMPLEX */


#ifndef __SUPERLU_SCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_SCOMPLEX

/* 
 * This header file is to be included in source files c*.c
 */
#ifndef SCOMPLEX_INCLUDE
#define SCOMPLEX_INCLUDE

typedef struct { float r, i; } lscomplex;


/* Macro definitions */

/*! \brief Complex Addition c = a + b */
#define c_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/*! \brief Complex Subtraction c = a - b */
#define c_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/*! \brief Complex-Double Multiplication */
#define cs_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/*! \brief Complex-Complex Multiplication */
#define cc_mult(c, a, b) { \
	float cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

#define cc_conj(a, b) { \
        (a)->r = (b)->r; \
        (a)->i = -((b)->i); \
    }

/*! \brief Complex equality testing */
#define c_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in scomplex.c */
void c_div(lscomplex *, lscomplex *, lscomplex *);
double c_abs(lscomplex *);     /* exact */
double c_abs1(lscomplex *);    /* approximate */
void c_exp(lscomplex *, lscomplex *);
void r_cnjg(lscomplex *, lscomplex *);
double r_imag(lscomplex *);
lscomplex c_sgn(lscomplex *);
lscomplex c_sqrt(lscomplex *);


#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_SCOMPLEX */

