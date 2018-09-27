// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include "ellipsoids.h"
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <random.hpp>
#include <random/uniform_int.hpp>
#include "simplex_samplers.h"
#include "copulas.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix copula_hypEll (Rcpp::NumericVector hyplane, Rcpp::NumericMatrix E,
                                   unsigned int num_slices, unsigned int numpoints){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef copula_ellipsoid<Point> CopEll;

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = hyplane.size(), i, j;

    std::vector<NT> hyp1(dim, 0.0);
    std::vector<std::vector<NT> > Gin(dim, std::vector<NT>(dim));

    for (i=0; i<dim; i++){
        hyp1[i] = hyplane[i];
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