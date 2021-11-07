#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/ballintersectsimplex.h"
#include "random_walks/gcw_estimator.hpp"
#include "random_walks/uniform_great_cycle_walk.hpp"
#include "check_convergence_test.h"


//' @export
// [[Rcpp::export]]
Rcpp::List check_convergence(Rcpp::NumericMatrix A,
                             Rcpp::NumericVector b,
                             Rcpp::NumericMatrix V,
                             Rcpp::NumericVector x0,
                             Rcpp::NumericMatrix samples,
                             unsigned int nu,
                             double lb,
                             double ub,
                             bool last_round)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef UnitBallIntersectSimplex <NT, VT, MT> Body;

    unsigned int d = A.ncol();
    VT Vnorms;
    bool too_few = false, precheck = false;
    NT ratio, alpha = 0.1;
    RNGType rng(d);

    Body BS(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), Rcpp::as<MT>(V), Rcpp::as<VT>(x0));
    MT X = Rcpp::as<MT>(samples);

    std::pair< std::pair<bool,bool>, std::pair<NT, VT> > res = check_convergence_test<VT>(BS,
                                                                                      X,
                                                                                      too_few,
                                                                                      ratio,
                                                                                      nu,
                                                                                      last_round,
                                                                                      alpha,
                                                                                      lb,
                                                                                      ub,
                                                                                      rng);
    
    return Rcpp::List::create(Rcpp::Named("conv") = res.first.first, Rcpp::Named("too_few") = res.first.second,
                              Rcpp::Named("ratio") = res.second.first, Rcpp::Named("x") = res.second.second);

}