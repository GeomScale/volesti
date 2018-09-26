// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.



#ifndef RSAMPLE_POINTS_H
#define RSAMPLE_POINTS_H

Rcpp::NumericMatrix Rsample_points (Rcpp::NumericMatrix A, unsigned int walk_len, double e, Rcpp::NumericVector InnerVec,
                           bool CG, bool ball_walk, double delta, bool Vpoly, bool Zono, bool sam_simplex,
                           bool sam_can_simplex, bool sam_arb_simplex, bool sam_ball, bool sam_sphere,
                           unsigned int numpoints, double variance);

#endif
