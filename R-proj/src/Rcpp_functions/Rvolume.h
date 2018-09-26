// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RVOLUME_H
#define RVOLUME_H

#include <Rcpp.h>

double Rvolume (Rcpp::NumericMatrix A, unsigned int walk_len, double e,
                Rcpp::NumericVector InnerVec, bool CG, unsigned int win_len,
                unsigned int N, double C, double ratio, double frac,
                bool ball_walk, double delta, bool Vpoly, bool Zono,
                bool coord, bool rounding);

#endif
