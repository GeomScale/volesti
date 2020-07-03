// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "diagnostics/raftery.hpp"

//'  An internal Rccp function for the random rotation of a convex polytope
//'
//' @param P A convex polytope (H-, V-polytope or a zonotope).
//' @param T Optional. A rotation matrix.
//' @param seed Optional. A fixed seed for the random linear map generator.
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A matrix that describes the rotated polytope
// [[Rcpp::export]]

Rcpp::NumericMatrix raftery(Rcpp::NumericMatrix samples)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT q = 0.025, r = 0.005, s = 0.95;

    std::pair<MT,VT> res = perform_raftery<VT>(Rcpp::as<MT>(samples), q, r, s);

    std::cout<<"res1 = "<<res.first<<std::endl;
    std::cout<<"res2 = "<<res.second<<std::endl;

    return Rcpp::wrap(res.first);

}
