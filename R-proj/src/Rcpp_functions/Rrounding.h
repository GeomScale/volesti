// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RROUNDING_H
#define RROUNDING_H

Rcpp::NumericMatrix rounding (Rcpp::NumericMatrix A, unsigned int walk_len, bool coord,
                              bool ball_walk, double delta, bool Vpoly, bool Zono);

#endif
