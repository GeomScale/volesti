// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#include <Rcpp.h>
#include <RcppEigen.h>
#include "exact_vols.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double SliceOfSimplex(Rcpp::NumericVector a, double z0){

    unsigned int dim = a.size();
    if (dim < 2) {
        throw Rcpp::exception("Dimension has to be greater than 2");
    }
    std::vector<double> hyp = Rcpp::as<std::vector<double> >(a);

    return vol_Ali(hyp,-z0, dim);

}