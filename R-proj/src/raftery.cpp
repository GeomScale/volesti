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
//' @param samples The sampled points from a geometric random walk.
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A matrix that describes the rotated polytope
// [[Rcpp::export]]

Rcpp::NumericMatrix raftery(Rcpp::NumericMatrix samples,
                            Rcpp::Nullable<double> q = R_NilValue,
                            Rcpp::Nullable<double> r = R_NilValue,
                            Rcpp::Nullable<double> s = R_NilValue)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT _q = (!q.isNotNull()) ? 0.025 : Rcpp::as<NT>(q);
    NT _r = (!r.isNotNull()) ? 0.05 : Rcpp::as<NT>(r);
    NT _s = (!s.isNotNull()) ? 0.95 : Rcpp::as<NT>(s);

    MT runs = Rcpp::as<MT>(samples).transpose();
    std::pair<MT,VT> res = perform_raftery<VT>(runs, _q, _r, _s);

    //std::cout<<"res1 = "<<res.first<<std::endl;
    //std::cout<<"res2 = "<<res.second<<std::endl;

    MT results(samples.rows(), 6);
    results << res.first, res.second;

    return Rcpp::wrap(results);

}
