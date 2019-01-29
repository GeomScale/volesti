// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include "ellipsoids.h"
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include "simplex_samplers.h"
#include "copulas.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix copula2 (Rcpp::NumericVector h, Rcpp::NumericMatrix E,
                                   Rcpp::Nullable<unsigned int> numSlices, Rcpp::Nullable<unsigned int> N){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef copula_ellipsoid<Point> CopEll;

    unsigned int num_slices = 100, numpoints = 4000000;

    if (numSlices.isNotNull()) {
        num_slices = Rcpp::as<unsigned int>(numSlices);
    }

    if (N.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = h.size(), i, j;

    std::vector<NT> hyp1 = Rcpp::as<std::vector<NT> >(h);;
    std::vector<std::vector<NT> > Gin(dim, std::vector<NT>(dim));

    for (i=0; i<dim; i++){
        hyp1[i] = h[i];
        for(j=0; j<dim; j++){
            Gin[i][j]=E(i,j);
        }
    }
    CopEll Ell(Gin);
    StdCopula = hypfam_ellfam<Point, RNGType >(dim, numpoints, num_slices, hyp1, Ell);

    for(i=0; i<num_slices; i++) {
        for(j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}