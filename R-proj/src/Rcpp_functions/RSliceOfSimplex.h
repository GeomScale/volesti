// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RSLICEOFSIMPLEX_H
#define RSLICEOFSIMPLEX_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double SliceSimplex(Rcpp::NumericVector hyplane){

    unsigned int dim = hyplane.size() - 1;
    NT z0 = hyplane[dim];
    std::vector<NT> hyp(dim, 0.0);

    for ( i=0; i<dim; i++) {
        hyp[i] = hyplane1[i];
    }

    return vol_Ali(hyp,-z0, dim);

}

#endif
