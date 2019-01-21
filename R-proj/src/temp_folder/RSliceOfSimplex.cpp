// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include "exact_vols.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double SliceSimplex(Rcpp::NumericVector hyplane){

    unsigned int dim = hyplane.size() - 1;
    double z0 = hyplane[dim];
    std::vector<double> hyp(dim, 0.0);

    for (unsigned int i=0; i<dim; i++) {
        hyp[i] = hyplane[i];
    }

    return vol_Ali(hyp,-z0, dim);

}