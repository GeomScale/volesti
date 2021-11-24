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
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/ballintersectsimplex.h"
#include "random_walks/fischer_gcw_walk.hpp"
#include "random_walks/uniform_great_cycle_walk.hpp"
#include "volume/volume_fischer_annealing_fast.hpp"

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
Rcpp::List compute_fischer_annealing(Rcpp::NumericMatrix A,
                                    Rcpp::NumericVector b,
                                    Rcpp::NumericVector mu_,
                                    Rcpp::NumericVector x0,
                                    unsigned int W)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef UnitBallIntersectSimplex <NT, VT, MT> Body;

    unsigned int d = A.ncol();
    VT b2 = Rcpp::as<VT>(b) - Rcpp::as<MT>(A)*Rcpp::as<VT>(x0);

    MT V2;

    Body BS(d, Rcpp::as<MT>(A), b2, V2, VT::Zero(d));

    //Body BS2(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), Rcpp::as<MT>(V), Rcpp::as<VT>(x0));

    //Body BS3(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), VV, VT::Zero(d), Vnorms_shifted);

    RNGType rng(d);
    VT p = Rcpp::as<VT>(mu_), center = Rcpp::as<VT>(x0);
    VT mu = Rcpp::as<VT>(mu_) - center;
    p -= center;

    std::pair<VT, VT> res = compute_annealing_fischer<GCWalk, FischerGCWalk, MT, NT>(BS, p, mu, rng, W);

    return Rcpp::List::create(Rcpp::Named("variances") = Rcpp::wrap(res.first), Rcpp::Named("ratios") = Rcpp::wrap(res.second));   
}
