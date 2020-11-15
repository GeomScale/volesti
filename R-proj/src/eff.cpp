// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20014-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "diagnostics/effective_sample_size.hpp"

//' Gelman-Rubin and Brooks-Gelman Potential Scale Reduction Factor (PSRF) for each marginal
//'
//' @param samples A matrix that contans column-wise the sampled points from a geometric random walk.
//' @param method A string to reauest diagnostic: (i) \code{'normal'} for psrf of Gelman-Rubin and (ii) \code{'interval'} for psrf of Brooks-Gelman.
//'
//' @references \cite{Gelman, A. and Rubin, D. B.,
//' \dQuote{Inference from iterative simulation using multiple sequences,} \emph{Statistical Science,} 1992.}
//'
//' @references \cite{Brooks, S. and Gelman, A.,
//' \dQuote{General Methods for Monitoring Convergence of Iterative Simulations,} \emph{Journal of Computational and Graphical Statistics,} 1998.}
//'
//' @return A vector that contains the values of PSRF for each coordinate
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eff(Rcpp::NumericMatrix samples)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
                              
    return Rcpp::wrap(eff_univariate<NT, VT>(Rcpp::as<MT>(samples)));
}

// [[Rcpp::export]]
Rcpp::NumericVector ESS_fft(Rcpp::NumericMatrix samples)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
                              
    return Rcpp::wrap(ess_univariate_fft<NT, VT>(Rcpp::as<MT>(samples)));
}

