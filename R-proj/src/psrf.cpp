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
#include "diagnostics/psrf.hpp"

//' Gelman-Rubin Potential Scale Reduction Factor (PSRF)
//'
//' @param samples A matrix that contans column-wise the sampled points from a geometric random walk.
//'
//' @references \cite{Gelman, A. and Rubin, D. B.,
//' \dQuote{Inference from iterative simulation using multiple sequences,} \emph{Statistical Science,} 1992.}
//'
//' @return The value of PSRF diagnostic.
// [[Rcpp::export]]
double psrf(Rcpp::NumericMatrix samples)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;    

    return perform_psrf<NT, VT>(Rcpp::as<MT>(samples));
}
