// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RVOLEXACT_H
#define RVOLEXACT_H

template <typename FT>
FT factorial(FT n);


double Rvol_exact (Rcpp::NumericMatrix A, bool exact_zono, bool exact_cube,
              bool exact_simplex, bool exact_cross, int dim);


#endif
