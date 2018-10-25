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
Rcpp::NumericMatrix copula_hyps (Rcpp::NumericVector hyplane1, Rcpp::NumericVector hyplane2,
                                 unsigned int num_slices, unsigned int numpoints){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = hyplane1.size(), i, j;

    std::vector<NT> hyp1(dim, 0.0);
    std::vector<NT> hyp2(dim, 0.0);
    typename std::vector<NT>::iterator hit1, hit2;
    i = 0;
    hit1 = hyp1.begin();
    hit2 = hyp2.begin();
    for ( ;  hit1!=hyp1.end(); ++hit1, ++hit2, ++i) {
        *hit1 = hyplane1[i];
        *hit2 = hyplane2[i];
    }

    StdCopula = twoParHypFam<Point, RNGType >(dim, numpoints, num_slices, hyp1, hyp2);

    for(i=0; i<num_slices; i++) {
        for(j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}