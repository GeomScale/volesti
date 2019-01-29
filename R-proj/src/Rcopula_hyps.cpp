// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
//#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include "simplex_samplers.h"
#include "copulas.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix copula1 (Rcpp::NumericVector h1, Rcpp::NumericVector h2,
                             Rcpp::Nullable<unsigned int> numSlices, Rcpp::Nullable<unsigned int> N){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    unsigned int num_slices = 100, numpoints = 4000000;

    if (numSlices.isNotNull()) {
        num_slices = Rcpp::as<unsigned int>(numSlices);
    }

    if (N.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = h1.size(), i, j;

    std::vector<NT> hyp1 = Rcpp::as<std::vector<NT> >(h1);
    std::vector<NT> hyp2 = Rcpp::as<std::vector<NT> >(h2);

    StdCopula = twoParHypFam<Point, RNGType >(dim, numpoints, num_slices, hyp1, hyp2);

    for(i=0; i<num_slices; i++) {
        for(j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}