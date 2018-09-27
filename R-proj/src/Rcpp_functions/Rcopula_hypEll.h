// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RCOPULA_HYPELL_H
#define RCOPULA_HYPELL_H


Rcpp::NumericMatrix copula_hypEll (Rcpp::NumericVector hyplane, Rcpp::NumericMatrix E,
                           unsigned int num_slices, unsigned int numpoints);

#endif
