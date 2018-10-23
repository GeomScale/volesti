/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998      Ross Ihaka
 *  Copyright (C) 2000-2013 The R Core Team
 *  Copyright (C) 2003      The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *   #include <Rmath.h>
 *
 *   double pnorm5(double x, double mu, double sigma, int lower_tail,int log_p);
 *         {pnorm (..) is synonymous and preferred inside R}
 *
 *   void   pnorm_both(double x, double *cum, double *ccum,
 *                     int i_tail, int log_p);
 *
 *  DESCRIPTION
 *
 *      The main computation evaluates near-minimax approximations derived
 *      from those in "Rational Chebyshev approximations for the error
 *      function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
 *      transportable program uses rational functions that theoretically
 *      approximate the normal distribution function to at least 18
 *      significant decimal digits.  The accuracy achieved depends on the
 *      arithmetic system, the compiler, the intrinsic functions, and
 *      proper selection of the machine-dependent constants.
 *
 *  REFERENCE
 *
 *      Cody, W. D. (1993).
 *      ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
 *      Special Function Routines and Test Drivers".
 *      ACM Transactions on Mathematical Software. 19, 22-32.
 *
 *  EXTENSIONS
 *
 *  The "_both" , lower, upper, and log_p  variants were added by
 *  Martin Maechler, Jan.2000;
 *  as well as log1p() and similar improvements later on.
 *
 *  James M. Rath contributed bug report PR#699 and patches correcting SIXTEN
 *  and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line
 *  with the original Cody code.
 */

/* Modified by Mutsuo Saito 2017-01-16 */

#include <math.h>
#include <float.h>

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

inline static double pnorm5(double x,
                            double mu,
                            double sigma,
                            int lower_tail,
                            int log_p)
{
    double p, cp;

    p = (x - mu) / sigma;
    x = p;
    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return(lower_tail ? p : cp);
}

inline static double log_pnorm(double x)
{

    return pnorm5(x, 0.0, 1.0, 1, 1);
}

inline static double pnorm(double x)
{

    return pnorm5(x, 0.0, 1.0, 1, 0);
}

inline static double log_dnorm(double x)
{

    return -0.5 * log(2 * M_PI) - 0.5 * x * x;
}

inline static double dnorm(double x)
{

    return exp(log_dnorm(x));
}
