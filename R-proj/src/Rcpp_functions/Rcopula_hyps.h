// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RCOPULA_HYPS_H
#define RCOPULA_HYPS_H

Rcpp::NumericMatrix copula_hyps (Rcpp::NumericVector hyplane1, Rcpp::NumericVector hyplane2,
                           unsigned int num_slices, unsigned int numpoints);

#endif
