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
Rcpp::NumericVector return_first_inside(Rcpp::NumericMatrix A,
                                        Rcpp::NumericVector b,
                                        Rcpp::NumericMatrix V,
                                        Rcpp::NumericVector x0,
                                        Rcpp::NumericMatrix samples)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef UnitBallIntersectSimplex <NT, VT, MT> Body;

    unsigned int d = A.ncol();
    RNGType rng(d);

    Body BS(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), Rcpp::as<MT>(V), Rcpp::as<VT>(x0));
    MT X = Rcpp::as<MT>(samples);
    int NN = X.cols();
    VT x = VT::Zero(d);

    for (int i=0; i<NN; i++)
    {
        x = X.col(i);
        if (BS.is_in(x)==-1)
        {
            break;
        }
    }
    
    return Rcpp::wrap(x);
}

